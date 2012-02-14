#!/usr/bin/perl
use strict;
use Getopt::Long;
use Cwd;
use Time::HiRes qw(time gettimeofday tv_interval);


my ($input,$debug);
my $self = bless {};
$self->{starttime} = time();
print STDERR "[init] ",time()-$self->{starttime}," secs...\n" if ($debug);

# $self->{work_dir} = $ENV{'HOME'}.'/plinko/data',
$self->{work_dir} = '/net/isilon7/nobackup/ensembl/avilella/other/other/other/other/plinko';
$self->{assembler} = 'mira';
$self->{assembler_exe}{mira} = $ENV{'HOME'}.'/plinko/assembler/mira/bin/mira';
$self->{assembler_exe}{trinity} = $ENV{'HOME'}.'/plinko/assembler/trinity/latest/Trinity.pl';
$self->{assembler_exe}{newbler26} = '/homes/avilella/newbler/Newbler_v2.6/bin/runAssembly';

$self->{lsf_options} = undef;
GetOptions(
	   'i|input:s' => \$input,
           'd|debug:s' => \$debug,
           'work_dir:s' => \$self->{work_dir},
           'assembler:s' => \$self->{assembler},
           'trinity_exe:s' => \$self->{assembler_exe}{trinity},
           'lsf_options:s' => \$self->{lsf_options},
           'retry:s' => \$self->{retry},
          );



if ($input =~ /^([SE]R\w\d+)/) { $input = 'http://www.ebi.ac.uk/ena/data/view/reports/sra/fastq_files/internal/' . $input; };
if ($input =~ /^http\:\/\/www\.ebi\.ac\.uk\/ena\/data\/view\/reports\/sra\/fastq_files\/internal\/([ES]R\w\d+)/) {
  $self->{project}{id} = $1;
  $self->{project}{url} = $input;
  $self->create_work_dir;
  $self->download_text_file;
  print STDERR "[download] ",time()-$self->{starttime}," secs...\n" if ($debug);
  $self->preprocess_text_file;
  $self->prepare_input_for_assembler;
  print STDERR "[prepare_input] ",time()-$self->{starttime}," secs...\n" if ($debug);
  $self->execute_assembler;
  print STDERR "[execute_assembler] ",time()-$self->{starttime}," secs...\n" if ($debug);
}

sub create_work_dir {
  my $self = shift;
  $self->{work_dir} = $self->{work_dir} . "/plinko_" . $self->{project}{id};
  my $cmd = 'mkdir -p '. $self->{work_dir};
  unless (0 == system($cmd)) { die ("cannot create work_dir\n: $cmd\n $!");}
  chdir($self->{work_dir});

  return;
}

sub preprocess_text_file {
  my $self = shift;
  open TEXTFILE, $self->{project}{id} or die $!;
  while (<TEXTFILE>) {
    chomp $_;
    if ($_ =~ /No results found/) {
      print STDERR "No results found. File probably still not uploaded\n"; exit 0;
    }
    next if ($_ =~ /^Study/);
    $DB::single=$debug;1;#??
    my ($Study,$Sample,$Experiment,$Run,$Organism,$Instrument_Platform,
        $Instrument_Model,$Library_Name,$Library_Layout,$Library_Source,
        $Library_Selection,$Run_Read_Count,$Run_Base_Count,$File_Name,
        $File_Size,$md5,$Ftp) = split("\t", $_);
    $self->{project}{filenames}{$File_Name}{file} = $self->{work_dir} . "/" . $File_Name;
    $self->{project}{filenames}{$File_Name}{Instrument_Platform} = $Instrument_Platform;
    $self->{project}{platforms}{$Instrument_Platform}{count}++;

    $self->download_ftp($Ftp);
  }

  return;
}

sub download_ftp {
  my $self      = shift;
  my $Ftp = shift;
  my $cmd = "wget -c ". $Ftp;
  unless (0 == system($cmd)) { die ("cannot download Ftp $Ftp\n: $cmd\n $!");}

  return;
}

sub download_text_file {
  my $self = shift;
  $self->{project}{text_file} = $self->{work_dir} . "/" . $self->{project}{id} . ".text_file";
  my $id = $self->{project}{id};
  my $cmd = "rm -f $id && wget ". $self->{project}{url} ." -o " . $self->{project}{text_file};
  unless (0 == system($cmd)) { die ("cannot download project text_file\n: $cmd\n $!");}
  $DB::single=$debug;1;
  return;
}

sub prepare_input_for_assembler {
  my $self = shift;
  $self->prepare_input_for_assembler_mira    if ($self->{assembler} =~ /mira|newbler26/);
  $self->prepare_input_for_assembler_trinity if ($self->{assembler} =~ /trinity/);
  unless ($self->{ready_for_assembly}) {die "cannot ready project for assembly: $!";};
  return;
}

sub prepare_input_for_assembler_mira {
  my $self = shift;
  $self->cleanup_input_file_mira unless (0 < $self->{retry});
  foreach my $File_Name (keys %{$self->{project}{filenames}}) {
    my $file = $self->{project}{filenames}{$File_Name}{file};
    $self->uncompress_and_concatenate_mira($File_Name);
  }

  $self->{ready_for_assembly} = 1;
  return;
}

sub prepare_input_for_assembler_trinity {
  my $self = shift;
  $self->cleanup_input_file_trinity unless (0 < $self->{retry});
  foreach my $File_Name (keys %{$self->{project}{filenames}}) {
    my $file = $self->{project}{filenames}{$File_Name}{file};
    $self->uncompress_and_concatenate_trinity($File_Name);
  }

  $self->{ready_for_assembly} = 1;
  return;
}

sub define_input_file_mira {
  my $self = shift;
  my $File_Name = shift;

  my $mira_platform = lc($self->translate_platform_mira($self->{project}{filenames}{$File_Name}{Instrument_Platform}));
  $self->{project}{input_file_mira} = $self->{work_dir} . "/" . $self->{project}{id} . "_in." . $mira_platform . ".fastq";

  return;
}

sub define_input_file_trinity {
  my $self = shift;
  my $File_Name = shift;

  $self->{project}{input_file_trinity} = $self->{work_dir} . "/" . $self->{project}{id} . ".trinity.fastq";

  return;
}

sub cleanup_input_file_mira {
  my $self = shift;
  my $file_regexp = $self->{work_dir} . "/*.fastq";
  my $cmd = "rm -f " . $file_regexp;
  unless (0 == system($cmd)) { die ("cannot cleanup old input_file $file_regexp\n: $cmd\n $!");}
  
  return;
}

sub cleanup_input_file_trinity {
  my $self = shift;
  my $file_regexp = $self->{work_dir} . "/*.trinity.fastq";
  my $cmd = "rm -f " . $file_regexp;
  unless (0 == system($cmd)) { die ("cannot cleanup old input_file $file_regexp\n: $cmd\n $!");}
  
  return;
}

sub uncompress_and_concatenate_mira {
  my $self = shift;
  my $File_Name = shift;
  my $file = $self->{project}{filenames}{$File_Name}{file};
  $self->define_input_file_mira($File_Name);

  my $method = 'gunzip';
  my $cmd = "$method -c ". $file ." >> " . $self->{project}{input_file_mira};
  print STDERR "[gunzip input_file $file] ",time()-$self->{starttime}," secs...\n" if ($debug);
  next if ($self->{retry});
  unless (0 == system($cmd)) { die ("cannot $method $file\n: $cmd\n $!");}

  return;
}

sub uncompress_and_concatenate_trinity {
  my $self = shift;
  my $File_Name = shift;
  my $file = $self->{project}{filenames}{$File_Name}{file};
  $self->define_input_file_trinity($File_Name);

  my $method = 'gunzip';
  my $cmd = "$method -c ". $file ." >> " . $self->{project}{input_file_trinity};
  print STDERR "[gunzip input_file $file] ",time()-$self->{starttime}," secs...\n" if ($debug);
  next if ($self->{retry});
  unless (0 == system($cmd)) { die ("cannot $method $file\n: $cmd\n $!");}

  return;
}

sub execute_assembler {
  my $self = shift;
  
  if ($self->{lsf_options} =~ /^(\d+)$/) {
    print STDERR "# $1\n";
    $self->{memory} = $1; my $mem = $self->{memory};
    $self->{lsf_options} = "bsub -I -q research-rh6 -R \"rusage[mem=$mem]\" -M$mem ";
  }

  $self->cmd_assembler_mira if ($self->{assembler} =~ /mira/);
  $self->cmd_assembler_newbler26 if ($self->{assembler} =~ /newbler26/);
  $self->cmd_assembler_trinity if ($self->{assembler} =~ /trinity/);
  my $cmd;
  my $work_dir = $self->{work_dir};
  print "\n cd $work_dir && ";
  if (defined $self->{lsf_options}) {
    $self->{project}{assembler_cmd} = $self->{lsf_options} . ' "' . $self->{project}{assembler_cmd} . '"';
  }
  $cmd = $self->{project}{assembler_cmd};

  print " $cmd\n";
#  unless (0 == system($cmd)) { die ("cannot create work_dir\n: $cmd\n $!");}

  return;
}

sub cmd_assembler_mira {
  my $self = shift;
  my $mira_exe = $self->{assembler_exe}{mira};

  my $mira_project_id = $self->{project}{id};
  my $prevalent_platform = undef;
  foreach my $this_platform (sort {$b->{count} <=> $a->{count}} keys %{$self->{project}{platforms}}) {
    my $count = $self->{project}{platforms}{$this_platform}{count};
    $prevalent_platform = $this_platform unless (defined $prevalent_platform);
    print STDERR "input: $this_platform - $count\n" if ($debug);
  }
  my $mira_platform = $self->translate_platform_mira($prevalent_platform);
  my $platform_settings;
  # $platform_settings = uc($mira_platform) ."_SETTINGS";
  # $platform_settings .= " -LR:mxti=no" if ($platform_settings =~ /solexa/i);
  my $cmd = "$mira_exe --project=$mira_project_id -DI:trt=/tmp/ --job=denovo,est,accurate,$mira_platform $platform_settings";
  $self->{project}{assembler_cmd} = $cmd;
}

sub cmd_assembler_newbler26 {
  my $self = shift;
  my $newbler26_exe = $self->{assembler_exe}{newbler26};
  my $project_id = $self->{project}{id};

  my $inputfile = $self->{project}{input_file_mira} || die "undefined newbler26 inputfile:$!";
  my $cmd = "$newbler26_exe -cdna -o newbler26.$project_id $inputfile";
  $self->{project}{assembler_cmd} = $cmd;
}


sub cmd_assembler_trinity {
  my $self = shift;
  my $trinity_exe = $self->{assembler_exe}{trinity};

  my $trinity_project_id = $self->{project}{id};
  my $inputfile = $self->{project}{input_file_trinity} || die "undefined trinity inputfile:$!";
  my $heap = '64G'; if (defined ($self->{memory})) { $heap = $self->{memory}/1000 . 'G';}
  my $cmd = "$trinity_exe --seqType fq --single $inputfile --CPU 1 --bflyHeapSpace $heap";
  $self->{project}{assembler_cmd} = $cmd;
}

sub translate_platform_mira {
  my $self = shift;
  my $platform = shift;
  $platform = 'solexa' if ($platform =~ /ILLUMINA/i);
  $platform = '454' if ($platform =~ /454/i);

  return $platform;
}
