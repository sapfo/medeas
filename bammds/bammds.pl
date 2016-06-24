#!/usr/bin/perl -w

# Copyright (C) 2013,2014 Ole Tange, Mike DeGiorgio, Anna-Sapfo
# Malaspinas, Jose Victor Moreno-Mayar, Yong Wang and Free Software
# Foundation, Inc.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License,
# or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 NAME

bammds - Creates a multidimensional scaling (MDS) plot of populations
for genetic data. The input files are a reference panel with genotypes
and one or several bam files.

=head1 SYNOPSIS

B<bammds> [--legend I<legendfile>] [--dim1 I<d1>] [--dim2 I<d2>] [--no-summary]
[--mapquality I<qual>] [--ntquality I<qual>] [--output I<outputfile>] [--pca] [--mds]
[I<file.bam> ...] [I<ref.bed> | I<ref.ped> | I<ref.tped> | I<ref.txt>
| I<ref.vcf> >] [--tfam I<file.tfam>] [--downsample I<proportion>]

B<bammds> [-l I<legendfile>] [-x I<d1>] [-y I<d2>] [--nosum] [-m I<qual>] [-n
I<qual>] [-o I<outputfile>] [--pca] [--mds] [I<file.bam> ...] [I<ref.bed> | I<ref.ped>
| I<ref.tped> | I<ref.txt> | I<ref.vcf> >] [--tfam I<file.tfam>] [-z I<proportion>]


=cut

use strict;
use local::lib;
use File::Temp qw(tempfile tempdir);
use IPC::Open3;
use Pod::Usage;
use File::Basename;
use Getopt::Long;
use File::Which;
use Text::CSV;

Getopt::Long::Configure("bundling");

# bammds_intersect --legend my_legend.csv --mapquality 30 --ntquality 20 --damagedread --undamagednt *.bam *.tped *.mike 

$Global::version = 20140602;
$Global::progname = "bammds";
my @retval = GetOptions
    ("debug|D" => \$opt::debug,
     "verbose|v" => \$opt::verbose,
     "help|h" => \$opt::help,
     "damagedread|d" => \$opt::damagedread,
     "undamagednt|u" => \$opt::undamagednt,
     "mapquality|m=i" => \$opt::mapq,
     "ntquality|n=i" => \$opt::ntq,
     "legend|legendfile|l=s" => \$opt::legend,
     "tfam=s" => \$opt::tfam,
     "dim1|xvector|x=s" => \$opt::xvector,
     "dim2|yvector|y=s" => \$opt::yvector,
     "triallelic_before_sampling|tribefore" => \$opt::triallelic_before_sampling,
     "outputfile|output|o=s" => \$opt::output,
     "mds" => \$opt::mds,
     "pca" => \$opt::pca,
     "no-summary|nosummary|nosum|no-sum" => \$opt::nosum,
     "downsample|z=f" => \$opt::subsample,
     "version|V" => \$opt::version,
    ) or pod2usage(2);
if($opt::help) { pod2usage(1) }
if($opt::version) {
  version();
  exit(0);
}
if($#ARGV < 0) {
  error("file.bam, file.tped, or file.mike is required. Please type man $Global::progname\n");
  pod2usage(-sections => "NAME");
}
if(defined($opt::ntq) and $opt::ntq > 40) {
  exit_error(1,"--ntquality must be below 40\n");
}
if(not which("bammds_plot.r")) {
  # Try adding the dir of $0 to the path
  my $exe_dir = $0;
  $exe_dir =~ s:/[^/]+$::;
  $ENV{'PATH'} .= ":".$exe_dir;
}
if(which("buffer")) {
  $::buffer_exe = "buffer";
} elsif(which("mbuffer")) {
  $::buffer_exe = "mbuffer";
} else {
  $::buffer_exe = "cat";
}

version();

# Use the locally created tmp-dir unless the user has set $TMPDIR
$ENV{'TMPDIR'} ||= "tmp";

my ($mike_file, $tped_arg_file, $ped_file, $bed_file, $vcf_file, @bam_files) = parse_file_args(@ARGV);

my ($reference_file, $tped_file, $chrPosRef, $chrVecRef, $posVecRef, @header) = 
  read_headers($mike_file, $tped_arg_file, $ped_file, $bed_file, $vcf_file, @bam_files);
my $asd_file = generate_asd($chrPosRef, $chrVecRef, $posVecRef, \@header, $reference_file, $mike_file, $tped_file, @bam_files);
my $legend = generate_legend($opt::legend, $header[0], $header[1], $reference_file, @bam_files);
my $plot_file = $opt::output;
plot_asd($asd_file, $legend, $plot_file, $opt::xvector, $opt::yvector);

sub parse_file_args {
  # Parse file name arguments
  # Input:
  #   @argv = @ARGV
  # Return:
  #   $mike_file = filename or undef
  #   $tped_file = filename or undef
  #   $ped_file = filename or undef
  #   $bed_file = filename or undef
  #   $vcf_file = filename or undef
  #   @bam_files = list of bams
  my @argv = @_;

  my ($mike_file, $tped_file, $ped_file, $bed_file, $vcf_file, @bam_files);
  
  for (@argv) {
    if(/.mike$/ or /.txt$/) {
      if($mike_file) {
	exit_error(1,"You can only have one .mike or .txt file\n");
      }
      $mike_file = $_;
    } elsif(/.tped$/) {
      if($tped_file) {
	exit_error(1,"You can only have one .tped file\n");
      }
      $tped_file = $_;
    } elsif(/.ped$/) {
      if($ped_file) {
	exit_error(1,"You can only have one .ped file\n");
      }
      $ped_file = $_;
    } elsif(/.bed$/) {
      if($bed_file) {
	exit_error(1,"You can only have one .bed file\n");
      }
      $bed_file = $_;
    } elsif(/.vcf$/) {
      if($vcf_file) {
	exit_error(1,"You can only have one .vcf file\n");
      }
      $vcf_file = $_;
    } elsif(/.bam$/) {
      push @bam_files, $_;
    } else {
      die("File type not known: $_");
    }
  }
  return ($mike_file, $tped_file, $ped_file, $bed_file, $vcf_file, @bam_files);  
}


sub plot_asd {
  my ($asd_file, $legend, $plot_file, $xvector, $yvector) = @_;
  my $legend_filled = fill_legend($legend);
  if(not defined($plot_file)) {
    my $mapQualityCutoff = $opt::mapq;
    my $ntQualityCutoff = $opt::ntq;
    my $quality_string = ($mapQualityCutoff ? ".".$mapQualityCutoff : "") . ($ntQualityCutoff ? ".".$ntQualityCutoff : "");
    my $newer;
    ($newer, $plot_file) = cache_name($reference_file, $legend_filled, @bam_files, 
				      combined_name(@bam_files, $reference_file) .
				      $quality_string . ".pdf");
    $plot_file =~ s:tmp/::;
  }
  if(-s $asd_file < 200) {
    ::exit_error(1,"ASD-file ($asd_file) too short. Not plotting\n");
  }
  $xvector ||= 1;
  $yvector ||= 2;

  # bammds_plot.r tmp/Han_700000reads_hg19,HanPop.HanIndv_70000,Wollstein.hg19.chr21.chr22.one_allele.30.20.asd tmp/Han_700000reads_hg19,HanPop.HanIndv_70000,Wollstein.hg19.chr21.chr22.one_allele.legend.filled.csv tmp/Han_700000reads_hg19,HanPop.HanIndv_70000,Wollstein.hg19.chr21.chr22.one_allele.30.20.pdf
  command_in_path_or_exit("bammds_plot.r");
  my $mds = $opt::mds ? "T" : "F";
  my $pca = $opt::pca ? "T" : "F";
  my $sum = $opt::nosum ? "F" : "T";
  if(not $opt::pca and not $opt::mds) {
    $mds = "T";
  }
  my $tmp = $ENV{'TMPDIR'};
  # When installing local packages $TMPDIR must be /tmp (donno why)
  $ENV{'TMPDIR'} = undef;
  debug("bammds_plot.r $asd_file $legend $legend_filled $plot_file $xvector $yvector $mds $pca $sum\n");
  print `bammds_plot.r $asd_file $legend $legend_filled $plot_file $xvector $yvector $mds $pca $sum`;
  $ENV{'TMPDIR'} = $tmp;
}


{
  my $longest_path;

  sub combined_name {
    my @bam_tped = @_;
    # Input:
    #   @bam_tped = paths to input files (bams + ref)
    # Output:
    #   $combined_name = combined_name of bams+tped
    #     my_dir/mypop.myindv.something.bam myindv2.bam ref.tped -> mypop.myindv,myindv2,ref
    #     If the name is too long: Use a hash value:
    #     my_dir/mypop.myindv.something.bam myindv2.bam ref.tped -> mypop.myindv,myindv2,+1,8bHyd
    
    sub hash_name {
      # Input:
      #   $len = length of hash value (up to 10)
      #   @names = data to make hash of
      my $len = shift;
      my @names = @_;
      my @values=("A".."Z",0..9,"_");
      use Digest::MD5 qw(md5 md5_hex);
      my ($val64bit) = unpack('Q', md5("@names"));
      my $mod = "";
      while($val64bit > 1 and $len) {
	$len--;
	$mod .= $values[ $val64bit % ($#values+1) ];
	$val64bit = $val64bit/($#values+1);
      }
      return $mod;
    }
    
    sub too_long_path {
      # Test if a file of length $i is too long
      # Input:
      #   $i = how long path to test
      my $i = shift;
      $ENV{"LANG"} = "C";
      open(my $dummy,"<","x"x$i);
      # We assume that we get an error string containing 'too long' if
      # the file name is too long.
      my $retval = ($! =~ /too long/);
      close $dummy;
      return $retval;
    }
    
    sub longest_path_sub {
      # Binary search for the longest path
      my ($lower,$upper) = @_;
      if($lower == $upper or $lower == $upper-1) { return $lower; }
      my $middle = int (($upper-$lower)/2 + $lower);
      if (too_long_path($middle)) {
	return longest_path_sub($lower,$middle);
      } else {
	return longest_path_sub($middle,$upper);
      }
    }
    
    sub longest_path {
      if(not $longest_path) {
	$longest_path = longest_path_sub(1,10000);
	debug("Longest filename: $longest_path\n");
      }
      return $longest_path;
    }
    
    sub basename_list {
      my @path = @_;
      for my $s (@path) {
	$s =~ s:^.*/([^/]+)/?$:$1:; # Remove dir from argument. If ending in /, remove final /
	$s =~ s:\.(bam|tped|ped|txt|text|bed)$::; # Remove .bam, .tped, .ped, .bed, .txt
      }
      return @path;
    }
    
    my @base = basename_list(@bam_tped);
    my $combined = "";
    my $j;
    for my $i (0 .. $#base) {
      # pop.indv.something.bam -> pop.indv
      my $n = $base[$i];
      $n =~ s/^([^.]+.[^.]+).*$/$1/;
      # length("+100,HASHED.legend.csv") = 22
      if(length($combined) + length($n) < longest_path() - 23) {
	$combined .= $n.",";
      } else {
	$j = $i;
	last;
      }
    }
    $combined =~ s/,$//;
    
    if(defined($opt::subsample)){
    	if($#bam_files==0){
			my $DSprop = $opt::subsample;
			$combined= "DS_" . $DSprop . "_" . $combined;
    	}else{
    		die("Error: Downsampling is only implemented for one sample")
    	}
    }

    if($j) {
      # +7,hashed
      my $missing = $#base - $j;
      my $hashname = hash_name(6,@base[$j .. $#base]);
      $combined .= "+$missing,$hashname";
    }
    return $combined;
  }
}

sub __READ_HEADERS__ {}


sub read_headers {
  # Input:
  #   $mike_file = filename or empty
  #   $tped_file = filename or empty
  #   $ped_file = filename or empty
  #   $bed_file = filename or empty
  #   $vcf_file = filename or empty
  #   @bam_files = list of bams
  # Return:
  #   $reference_file = filename that was not empty
  #   $tped_file = filename of generated .tped if any
  #   $chrPosRef = \%chrPos = $chrPos{$chromo} is list of positions for this chromosome
  #   $chrVecRef = \@chrVec = list of chromosomes
  #   $posVecRef = \@posVec = list of positions (matching @chrVec)
  #   @header = ("Individual1 \t individual2", "population_of_indv1 \t population_of_indv2")

  my ($mike_file, $tped_file, $ped_file, $bed_file, $vcf_file) = @_;
  my ($reference_file, $chrPosRef, $chrVecRef, $posVecRef, @header);
  if($mike_file) {
    $reference_file = $mike_file;
    ($chrPosRef,$chrVecRef,$posVecRef) = mikeformat_to_chrpos($mike_file);
    if($opt::tfam) {
      @header = tped_population_individual_header($opt::tfam);
    } else {
      @header = mikeformat_population_individual_header($mike_file)
    }
  } elsif($tped_file) {
    $reference_file = $tped_file;
    ($chrPosRef,$chrVecRef,$posVecRef) = tped_to_chrpos($tped_file);
    @header = tped_population_individual_header($tped_file);
  } elsif($ped_file) {
    # Convert ped to tped and run the rest of the script as if you had a .tped
    $reference_file = $ped_file;
    $tped_file = ped_to_tped($ped_file);
    ($chrPosRef,$chrVecRef,$posVecRef) = tped_to_chrpos($tped_file);
    @header = tped_population_individual_header($tped_file);
  } elsif($bed_file) {
    # Convert bed to tped and run the rest of the script as if you had a .tped
    $reference_file = $bed_file;
    $tped_file = bed_to_tped($bed_file);
    ($chrPosRef,$chrVecRef,$posVecRef) = tped_to_chrpos($tped_file);
    @header = tped_population_individual_header($tped_file);
  } elsif($vcf_file) {
    # Convert vcf to tped and run the rest of the script as if you had a .tped
    $reference_file = $vcf_file;
    $tped_file = vcf_to_tped($vcf_file);
    ($chrPosRef,$chrVecRef,$posVecRef) = tped_to_chrpos($tped_file);
    @header = tped_population_individual_header($tped_file);
  } else {
    exit_error(1,"Genotypes must be delivered in .ped, .bed, .tped, .vcf or .mike format\n");
  }

  # Header from bam files
  for my $file (@bam_files) {
    my ($pop,$indv);
    my $bam = $file;
    # foo/bar.bam => bar.bam
    $bam =~ s{.*/}{}g;
    if($bam =~ m/([^\.]+)\.([^\.]+)(\..*)?\.bam$/i) {
      # population.individual.garbage.bam
      # population.individual.bam
      $pop = $1;
      $indv = $2;
    } elsif($bam =~ /^([^\.]+)\.bam$/i) {
      # individual.bam
      $pop = $1;
      $indv = $1;
    } else {
      die_bug("Cannot compute individual from bam name");
    }
    # Prepend the individual and population to the header
    $header[0] = $indv."\t".$header[0];
    $header[1] = $pop."\t".$header[1];
  }

  return ($reference_file, $tped_file, $chrPosRef, $chrVecRef, $posVecRef, @header);
}

sub mikeformat_population_individual_header {
  # Input:
  #   $mike_file = filename of .mike
  # Output:
  #   "indv1 \t indv2","pop1 \t pop2"
  my $mike_file = shift;
  open(IN,"<",$mike_file) || exit_error(1,"Cannot open $mike_file\n");
  my $indv = <IN>;
  chomp $indv;
  my $pop = <IN>;
  chomp $pop;
  close IN;
  return ($indv,$pop);
}


sub tped_population_individual_header {
  # Input:
  #   $tped_file = filename of tped
  # Output:
  #   "indv1 \t indv2","pop1 \t pop2"
  my $tped_file = shift;
  if($opt::tfam) {
    $tped_file = $opt::tfam;
  } else {
    $tped_file =~ s/\.tped$/\.tfam/ or exit_error(1,".tped file must have a .tfam partner\n");
  }
  open(IN,"<",$tped_file) || exit_error(1,".tped file must have a .tfam partner\n");
  my @indv_header;
  my @pop_header;
  my $last_pop = "";
  while(<IN>) {
    # "POL"   "Pl26"  0       0       0       0.7
    my ($pop,$indv) = split/\s+/,$_,3;
    if($pop eq "" or $pop eq '""') {
      # Empty means last
      $pop = $last_pop;
    } else {
      $last_pop = $pop;
    }
    push @indv_header, $indv;
    push @pop_header, $pop;
  }
  close IN;
  return ((join "\t", @indv_header),(join "\t", @pop_header));
}


sub tped_to_chrpos {
  my $tped_file = shift;
  open(IN, "<", $tped_file) || exit_error(1,"Cannot open $tped_file\n");
  my @chrVec = ();
  my @posVec = ();
  my $chrPos;
  while(<IN>){
    # chr,pos_name, ?, pos, indv-1-parent-1, indv-1-parent-2, indv-2-parent-1, indv-2-parent-2, ...
    # 21 rs12627229 0 10913441 T T T T T T T T T T T T T T C C T T 0 0
    # ignore column 2,3,5
    my ($chr,$i2,$i3,$pos,$i5) = split /\s+/, $_, 5;
    # Map chrX => X, chromX => X
    $chr =~ s/^chr(om)?//i;
    push @{$chrPos->{$chr}},$pos;
    push @chrVec, $chr;
    push @posVec, $pos;
  }
  close(IN);
  return ($chrPos,\@chrVec,\@posVec);
}


sub mikeformat_to_chrpos {
  my $mike_file = shift;

  open(IN, "<", $mike_file) || exit_error(1,"Cannot open $mike_file\n");
  my $line = <IN>; # Individual IDs
  $line = <IN>; # Population IDs
  my @chrVec = ();
  my @posVec = ();
  my $chrPos;
  while(<IN>){
    # Mikes format:
    #   21 rs12627229 10913441 T C 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1
    #   chr21 rs12627229 10913441 T C 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1
    if(not /^(\w+)\s+\S+\s+(\d+)\s+/) {
      die;
    }
    my($chr,$pos) = ($1,$2);
    # Map chrX => X, chromX => X
    $chr =~ s/^chr(om)?//i;
    push @{$chrPos->{$chr}},$pos;
    push @chrVec, $chr;
    push @posVec, $pos;
  }
  close(IN);
  return ($chrPos,\@chrVec,\@posVec);
}


sub ped_to_tped {
  my $ped_file = shift;
  my $no_extension_ped = $ped_file;
  $no_extension_ped =~ s/.ped$//i;
  my ($newer,$tped_file) = cache_name($ped_file,$no_extension_ped.".tped");
  $tped_file =~ s/.tped$//i;
  if($newer) {
    # The .ped is newer => convert again
    shell_run_or_die("p-link","--file $no_extension_ped --out $tped_file --recode --transpose 2>&1");
  }
  return $tped_file.".tped";
}


sub bed_to_tped {
  my $bed_file = shift;
  my $no_extension_bed = $bed_file;
  $no_extension_bed =~ s/.bed$//i;
  my ($newer,$tped_file) = cache_name($bed_file,$no_extension_bed.".tped");
  $tped_file =~ s/.tped$//i;
  if($newer) {
    # The .bed is newer => convert again
    shell_run_or_die("p-link", "--bfile $no_extension_bed --out $tped_file --recode --transpose 2>&1");
  }
  return $tped_file.".tped";
}


sub vcf_to_tped {
  my $vcf_file = shift;
  my $no_extension_vcf = $vcf_file;
  $no_extension_vcf =~ s/.vcf$//i;
  my ($newer,$tped_file) = cache_name($vcf_file,$no_extension_vcf.".tped");
  $tped_file =~ s/.tped$//i;
  if($newer) {
    shell_run_or_die("vcftools", "--vcf $vcf_file --plink-tped --out $tped_file 2>&1");
  }
  return $tped_file.".tped";
}


sub __GENERATE_ASD__ {}


sub generate_asd {
  # If source files are newer or no asd file: Generate asd file
  # Input:
  #   $chrPosRef = \%chrPos = $chrPos{$chromo} is list of positions for this chromosome
  #   $chrVecRef = \@chrVec = list of chromosomes
  #   $posVecRef = \@posVec = list of positions (matching @chrVec)
  #   $header_ref = \@header = ("Individual1 \t individual2", "population_of_indv1 \t population_of_indv2")
  #   $reference_file = .mike/.tped/.ped/.bed file name
  #   $mike_file = reference file in mike format if any
  #   $tped_file = generated tped-file if any
  #   @bam_files = list of bam filenames
  # Return:
  #   $asd_file = name of asd file
  my ($chrPosRef, $chrVecRef, $posVecRef, $header_ref, $reference_file, $mike_file, $tped_file, @bam_files) = @_;
  my (@header) = @$header_ref;
  my $mapQualityCutoff = $opt::mapq;
  my $ntQualityCutoff = $opt::ntq;
  my $quality_string = ($mapQualityCutoff ? ".".$mapQualityCutoff : "") . ($ntQualityCutoff ? ".".$ntQualityCutoff : "");
  my ($newer, $asd_file) = cache_name($reference_file, @bam_files, 
				      "tmp/".combined_name(@bam_files, $reference_file) .
				      $quality_string . ".asd");
  my $chromosome_prefix = "";
  if($newer or -s $asd_file < 200) {
    # The source files are newer or asd_file is empty: Generate the asd-file.
    # Test if .asd is writable
    open(TEST_WRITABLE, ">", $asd_file) || die("Cannot open $asd_file");
    close TEST_WRITABLE;

    my ($sitesFile, $empty_bam);
    if(@bam_files) {
      # There is at least 1 bam file
      # Generate sites_file and dummy_bam-file
      ($chromosome_prefix,$empty_bam) = dummy_bam($bam_files[0],$chrVecRef,$posVecRef);
      $sitesFile = sites_file($chromosome_prefix,$chrVecRef,$posVecRef);
    }

    # Pipe to parallel merge
    command_in_path_or_exit(qw(bammds_parallel bammds_merge bammds_allelesharing bammds_allelesharingsum));
    my $debug_tee = $opt::debug ? "tee ${asd_file}_tee |" : "";
    my $tribefore = $opt::triallelic_before_sampling ? "--tribefore" : "";
    my $debug = $opt::debug ? "--debug" : "";
    my $verbose = $opt::verbose ? "--verbose" : "";
    my $parallel_cmd = qq(bammds_parallel -j70% --compress --round-robin --header 2 --block 500k --pipe '$::buffer_exe | bammds_merge $debug $verbose $tribefore | bammds_allelesharing');
    my $cmd = "$debug_tee $parallel_cmd | bammds_allelesharingsum >$asd_file";
    open(ASD, "|-", $cmd) || die("Cannot open $cmd");
    debug($cmd,"\n");
    # Print headers.
    print ASD map { $_,"\n" } @header;

    # Call samtools
    my $flagUseDamageRead = $opt::damagedread;
    my $flagUseUndamagedNT = $opt::undamagednt;

    if(@bam_files) {
      $mapQualityCutoff ||= 20;
      $ntQualityCutoff ||= 30;
	if(defined($opt::subsample)){
		my $DSprop = $opt::subsample;
		if($#bam_files==0){
			my $downsamplemessage="WARNING: Downsampled " . $DSprop . " from the original BAM file. This feature has only been tested with samtools 0.1.18 (r982:295). Some problems have been reported when proportions greater than 0.5 are used. Make sure that samtools view -bs works for you. \n";
			print STDERR $downsamplemessage;
			open(PU,"samtools view -bs $DSprop @bam_files | samtools mpileup -Bl $sitesFile -q $mapQualityCutoff -Q $ntQualityCutoff - $empty_bam | $::buffer_exe |") || die;
			debug("samtools view -bs $DSprop @bam_files | samtools mpileup -Bl $sitesFile -q $mapQualityCutoff -Q $ntQualityCutoff - $empty_bam | $::buffer_exe |\n");
		}else{
			die("Error: Downsampling is only implemented for one sample")
		}
	}else{
      open(PU,"samtools mpileup -Bl $sitesFile -q $mapQualityCutoff -Q $ntQualityCutoff @bam_files $empty_bam | $::buffer_exe |") || die;
      debug("samtools mpileup -Bl $sitesFile -q $mapQualityCutoff -Q $ntQualityCutoff @bam_files $empty_bam | $::buffer_exe |\n");
	}
    }
    my @tped_fh = open_tped_files($tped_file);
    if($mike_file) {
      open(MIKE, "<", $mike_file) || exit_error(1,"Cannot open $mike_file\n");
      my ($dummy);
      # Read 2 headers
      $dummy = <MIKE>;
      $dummy = <MIKE>;
    }
    my $i;
    for($i = 0; $i < $#$chrVecRef+1; $i++) {
      my @out;
      push @out, $chrVecRef->[$i], $posVecRef->[$i];
      if(@bam_files) {
	if(eof(PU)) { 
	  unlink $asd_file;
	  close PU;
	  close ASD;
	  error("The mpileup command does not work. ".
		"Are the chromosomes named the same way in all files (e.g. NN and not chrNN)?\n".
		"Are the chromosomes ordered the same way in all the files?\n");
	  for my $bam (@bam_files) {
	     error("BAM header of $bam:\n" . qx{ samtools view -H $bam });
	  }
	  exit_error(1,"Header of panel:\n" . qx{ samtools view -H $empty_bam });
	}
	$_ = <PU>;
	chomp;
	push @out, "PU", $_;
      }
      if(@tped_fh) {
	# Paste the tped files here.
	for my $fh (@tped_fh) {
	  # read a line
	  my $tped_line = <$fh>;
	  chomp $tped_line;
	  push @out, "TPED", $chromosome_prefix.$tped_line;
	}
      }
      if($mike_file) {
	my $mike_line = <MIKE>;
	chomp($mike_line);
	push @out, "Mike", $mike_line;
      }
      print ASD join("\0", @out),"\n";
    }
    close(PU);
    map { close $_ } @tped_fh;
    $opt::debug or $empty_bam and unlink $empty_bam;
    $opt::debug or $sitesFile and unlink $sitesFile;
    ###Fill in the last markers that are expected to be 0
    my $nMarkers = $#$chrVecRef+1;
    debug("Is $i < $nMarkers\n");
    if($i < $nMarkers) {
      print STDERR "Filling missing markers\n";
      die_bug("All markers should be in all files");
      my $nInds = $#bam_files+1;
      while($i < $nMarkers){
	my @out;
	push @out, $chrVecRef->[$i], $posVecRef->[$i];
	for(my $individual=1; $individual <= $nInds; $individual++) {
	  push @out, (0,0,0,0);
	}
	print ASD join("\t", @out),"\n";
	$i++;
      }
    }
    debug("close ASD\n");
    close(ASD);
    debug("closed ASD\n");
  } else {
    warning("Using pre-computed $asd_file\n");
    if(not $opt::nosum) {
	warning("Consider using --no-summary to speed up plotting and use less memory.\n");
    }
  }
  debug("Generated $asd_file\n");

  return $asd_file;
}


sub open_tped_files {
  my @filenames = @_;
  my @tped_fh;
  for (@filenames) {
    if(defined $_) {
      open(my $fh, "<", $_) || die;
      push @tped_fh, $fh;
    }
  }
  return @tped_fh;
}


sub max {
    # Returns:
    #   Maximum value of array
    my $max;
    for (@_) {
        # Skip undefs
        defined $_ or next;
        defined $max or do { $max = $_; next; }; # Set $_ to the first non-undef
        $max = ($max > $_) ? $max : $_;
    }
    return $max;
}


sub bam_header {
  my $inputbamfile = shift;
  command_in_path_or_exit("samtools");
  return `samtools view -H $inputbamfile`;
}


sub dummy_bam {
  # Generate a bam file with all the positions, so that samtools will
  # spit out all positions.
  #
  # Input:
  #   $inputbamfile = filename of existing BAM to grab the header from
  #   \@chr = list(1,1,1,2,2,2,2) - list of chromosomes matching the position in @pos
  #   \@pos = list(123000,133000,133300,100,23000,499949,600000)
  # Output:
  #   $chromosome_prefix = the prefix string for chromosomes in the bam (e.g. chr or chrom)
  #   $dummy_bam_file = bam file name with all chr,pos in it.
  my $inputbamfile = shift;
  my $chr_ref = shift;
  my $pos_ref = shift;
  my @chr = @$chr_ref;
  my @pos = @$pos_ref;
  my @dummy_sam;
  my @bam_header = bam_header($inputbamfile);
  push @dummy_sam, @bam_header;
  my $chromosome_prefix = chromosome_prefix_from_bam_header(@bam_header);
  my $readname = "r";
  my $samflag = 0;
  my $mapq = 255; # Maximal mapping quality: We donot want to lose this due to samtools -q limit.
  my $cigar = "1M";
  my $mate = "*";
  my $matepos = 0;
  my $templatelen = 0;
  my $readseq = "A";
  # Give it the maximal readquality. Q40
  my $readqual = "J";
  # Check if the input is sorted
  my $last_chr_pos = 0;
  my $chr_pos = 0;
  for (my $i = 0; $i <= $#chr; $i++) {
    my($chr,$pos) = ($chr[$i],$pos[$i]);
    $chr =~ s/chr(om)?//i; # Remove chr if it is there.
    $chr_pos = $chr*200_000_000_000 + $pos; # Express position as a single integer
    if($chr_pos < $last_chr_pos) {
	exit_error(1,"Your reference panel must be sorted by chromosome and position. It is not at chr $chr, pos $pos.\n",
		   "You may be able to fix it by doing something like:\n",
		   "  cat reference.txt | head -n 2 > ref.out\n",
		   "  cat reference.txt | tail -n +3 | sort -k1,1n -k3,3n >> ref.out\n",
		   "  mv ref.out reference.txt\n");
    }
    $last_chr_pos = $chr_pos;
    $chr = $chromosome_prefix.$chr; # Add the prefix detected from BAM-file
    push @dummy_sam, join("\t",$readname, $samflag, $chr, $pos, $mapq, $cigar,
			  $mate, $matepos, $templatelen, $readseq, $readqual),"\n";
  }
  if($opt::debug) {
    open (OUT, ">tmp/test.sam") || die("cannot write tmp/test.sam");
    print OUT @dummy_sam;
    close OUT;
  }
  my ($fh_empty, $dummy_bam_file) = ::tempfile(SUFFIX => ".bam");
  print $fh_empty sam_to_bam(@dummy_sam);
  close $fh_empty;
  return $chromosome_prefix,$dummy_bam_file;
}


sub chromosome_prefix_from_bam_header {
  # Input:
  #   @bam_header = header of BAM files
  # Output:
  #   $prefix = the prefix for each chromosome (e.g. chr, chrom, Chrom)
  my @bam_header = @_;
  my @chr;
  my $max_len = 0;
  my $prefix = "";
  for (@bam_header) {
    if(/SQ\s+SN:(\S+)/) {
      # @SQ     SN:chr1 LN:247249719
      # @SQ     SN:chrX LN:154913754
      my $chr = $1;
      # The SN can also be other stuff such as:
      # @SQ     SN:GL000229.1   LN:19913
      $chr =~ /chr|^c.\d\d/i or next;
      push @chr, $1;
      $max_len = max(length $1,$max_len);
    }
  }
  # Determine prefix by frequency goes wrong on SN:GL000192.1 .. SN:GL000229.1
  for(my $pos = 0; $pos < $max_len; $pos++) {
    my %freq;
    map { $freq{substr($_,$pos,1)}++ } @chr;
    my @top = sort { $freq{$b} <=> $freq{$a} } keys %freq;
    if(not defined $top[0]
       or
       $top[0] =~ /0-9/
       or
       $freq{$top[0]} < 0.50 * $#chr
      ) {
      # The most frequent is undef or a number or fewer than 50% => We are done
      last;
    }
    $prefix .= $top[0];
  }
  return $prefix;
}


sub sam_to_bam {
  # Input:
  #   @sam = content of samfile
  # Output:
  #   $bamfile = content of bamfile
  my @sam = @_;
  my $bam = "";
  my($wtr, $rdr, $err);
  use Symbol 'gensym'; 
  $err = gensym;
  my $pid = open3($wtr, $rdr, $err, "samtools view -S -b - 2>/dev/null");
  if(fork()) {
    # Parent
    close $wtr;
    local $/ = undef;
    while(<$rdr>) {
      $bam .= $_;
    }
    while(<$err>) {
      print $_;
    }
  } else {
    close $rdr;
    close $err;
    print $wtr @sam;
    exit(0);
  }
  return $bam;
}


sub sites_file {
  my $chromosome_prefix = shift;
  my $chrVecRef = shift;
  my $posVecRef = shift;
  my @chrVec = @$chrVecRef;
  my @posVec = @$posVecRef;
  my @out;

  for (my $i = 0; $i <= $#chrVec; $i++) {
    push @out, $chromosome_prefix.$chrVec[$i], "\t", $posVec[$i]-1, "\t", $posVec[$i], "\n";
  }

  my ($fh_sites, $sites_file) = ::tempfile(SUFFIX => ".bed");
  print $fh_sites @out;
  close $fh_sites;
  return $sites_file
}


sub __LEGEND__ {}


sub generate_legend {
  my ($legend_file, $individuals_tab, $populations_tab, $reference_file, @bam_files) = @_;
  my (@sample_populations, @sample_individuals,
      @reference_populations, @reference_individuals, @individuals,
      @populations);

  @populations = map { s/^"(.*)"/$1/; $_ } split /\s/, $populations_tab;
  @individuals = map { s/^"(.*)"/$1/; $_ } split /\s/, $individuals_tab;

  @sample_populations = @populations[0..$#bam_files];
  @sample_individuals = @individuals[0..$#bam_files];

  @reference_populations = @populations[$#bam_files+1..$#populations];
  @reference_individuals = @individuals[$#bam_files+1..$#individuals];

  my @table;
  # Header
  push @table, [ qw(Population pop_label Individual indv_label sample order Color pch cex Other columns are ignored) ],
    [ qw(* * * * * * * * *) ];

  # Populations - sample first
  my $first_sample = 1;
  my %seen;
  for my $p (@sample_populations) {
    $seen{$p}++ and next;
    push @table, [ $p,"","","",$first_sample ];
    $first_sample = "";
  }
  my $first_reference = 0;
  for my $p (@reference_populations) {
    $seen{$p}++ and next;
    push @table, [ $p,"","","",$first_reference ];
    $first_reference = "";
  }
  push @table, [ "" ];

  # Individuals
  my $first_individual = "*";
  for (my $i = 0; $i <= $#individuals; $i++) {
    push @table, [ $populations[$i], "NA", 
		   $individuals[$i], $first_individual, 
		   $first_individual, $first_individual,
		   $first_individual, $first_individual,
		   $first_individual ];
    $first_individual = "";
  }

  my @content = <<'_EOS' ;
# Lines starting with # are ignored

# Population: Internal ID for population. * = all populations
# pop_label: Label that will be used for the population in plots.
# pop_label: Population: * = use population as pop_label.
# pop_label: Individual: * = pop_label of Population.
# Individual: Internal ID for individual. * = all individuals
# indv_label: Label that will be used for the individual in plots.
# indv_label: * = Use individual as indv_label.
# sample: 1=sample (bam). 0=reference (tped or similar). -1=remove from plot.
# sample: *=use value defined for population
# order: Sort order of individuals in in ASD-file. No duplicates allowed.
# Color: Color given as hex-value.
# Color: See http://www.w3schools.com/tags/ref_colorpicker.asp
# Color: * means auto select based on pop_label/indv_label
# pch: R's pch-value. * means auto select
# pch: square brackets around a value means plot circle around the value
# cex: R's cex-value. For population * means auto select
# cex: For individual: Use value for population

# A specified individual/population takes priority over a *
# An empty field means the same as the field above
_EOS

  push @content, map { join(",", @$_) . "\n" } @table;

  mkdir "tmp";
  if(not $legend_file) {
    $legend_file = "tmp/".combined_name(@bam_files, $reference_file).".legend.csv";
  }
  if(not -e $legend_file or -s $legend_file < 30) {
    # Do not overwrite legend file
    my $out;
    open($out, ">", $legend_file) || exit_error(1,"Cannot write $legend_file: $!\n");
    print $out @content;
    close $out;
  } else {
    warning("Not overwriting $legend_file\n");
  }

  return $legend_file;
}


sub fill_legend {
  my ($legend) = @_;

  my $legend_file_no_csv = $legend;
  $legend_file_no_csv =~ s/.csv$//;
  my ($newer,$legend_filled) = cache_name($legend,$legend_file_no_csv.".filled.csv");
  if($newer or -s $legend_filled < 30) {
    debug("bammds_fillout $legend > $legend_filled\n");
    command_in_path_or_exit("bammds_fillout");
    print `bammds_fillout $legend > $legend_filled`;
  }
  return $legend_filled;
}


sub __GENERIC_LIBRARIES__ {}


sub my_dump {
    # Returns:
    #   ascii expression of object if Data::Dump(er) is installed
    #   error code otherwise
    my @dump_this = (@_);
    eval "use Data::Dump qw(dump);";
    if ($@) {
        # Data::Dump not installed
        eval "use Data::Dumper;";
        if ($@) {
            my $err =  "Neither Data::Dump nor Data::Dumper is installed\n".
                "Not dumping output\n";
            print $Global::original_stderr $err;
            return $err;
        } else {
            return Dumper(@dump_this);
        }
    } else {
	# Create a dummy Data::Dump:dump as Hans Schou sometimes has
	# it undefined
	eval "sub Data::Dump:dump {}";
        eval "use Data::Dump qw(dump);";
        return (Data::Dump::dump(@dump_this));
    }
}


sub shell_run_or_die {
  # Run the shell command.
  # If it fails: Print the output on STDERR and exit.
  #
  # Input:
  #   @cmd = shell command
  # Outout:
  #   N/A
  my @cmd = @_;
  debug(@cmd,"\n");
  command_in_path_or_exit($cmd[0]);
  my @out = `@cmd`;
  if($? or grep /error/i, @out) {
    error("This failed: @cmd\n");
    print STDERR @out;
    exit(1);
  }
}


sub command_in_path_or_exit {
  my @cmds = @_;
  for(@cmds) {
    which($_) or exit_error(1,"Command not found: $_\n");
  }
}


sub cache_name {
  # Given n paths:
  # return a name to a file in tmp/ based on $dst and whether the tmp-file is is newer than the first n-1
  #
  # Input:
  #   @src = full paths to source files
  #   $dst = full path to the would-be destination file
  # Output:
  #   $newer = bool: True if one of @src is newer than $tmp or $tmp does not exist
  #   $tmp = relative path to the cache file
  my (@src) = @_;
  my $dst = pop @src;
  mkdir ("tmp");
  my $tmp = "tmp/".basename($dst);
  my $newer = ! -e $tmp;
  if(not $newer) {
    for my $src (@src) {
      # -M = time since modification. 0 = now, 1.0 = 24 hours ago.
      $newer = $newer || (-M $tmp > -M $src);
    }
  }
  return ($newer,$tmp);
}


sub uniq {
  # Remove duplicates
  # Input:
  #   @values = list of values with duplicates
  # Return:
  #   @uniq = unique values
  return keys %{{ map { $_ => 1 } @_ }};
}

sub die_bug {
  my $bugid = shift;
  $Global::progname ||= $0;
  $Global::version ||= 0.1;
  print STDERR
    ("$Global::progname: This should not happen. You have found a bug.\n",
     "Please contact <parallel\@gnu.org> and include:\n",
     "* The version number: $Global::version\n",
     "* The bugid: $bugid\n",
     "* The command line being run\n",
     "* The files being read (put the files on a webserver if they are big)\n",
     "\n",
     "If you get the error on smaller/fewer files, please include those instead.\n");
  ::wait_and_exit(255);
}


sub exit_error {
  my $exit_value = shift;
  error(@_);
  exit $exit_value;
}


sub error {
    my @w = @_;
    my $fh = $Global::original_stderr || *STDERR;
    my $prog = $Global::progname || $0;
    print $fh $prog, ": Error: ", @w;
}


sub warning {
    my @w = @_;
    my $fh = $Global::original_stderr || *STDERR;
    my $prog = $Global::progname || $0;
    print $fh $prog, ": Warning: ", @w;
}


sub debug {
    # Returns: N/A
    $opt::debug or return;
    @_ = grep { defined $_ ? $_ : "" } @_;
    print STDERR @_;
}


sub version {
    # Returns: N/A
    print join("\n",
               "$Global::progname $Global::version",
               "Copyright (C) 2013,2014 Ole Tange, Mike DeGiorgio, Anna-Sapfo Malaspinas, Jose Victor Moreno-Mayar, Yong Wang and Free Software Foundation, Inc.",
               "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>",
               "This is free software: you are free to change and redistribute it.",
               "$Global::progname comes with no warranty.",
               "",
               "Web site: http://www.gnu.org/software/${Global::progname}\n",
	       "Send bugreports, beers and margaritas to: bammds-users\@nongnu.org\n",
	       "When using $Global::progname for a publication please cite:\n",
	       "Malaspinas A-S, Tange O, Moreno-Mayar JV, Rasmussen M, DeGiorgio M, Wang Y, et al. bammds: a tool for assessing the ancestry of low-depth whole-genome data using multidimensional scaling (MDS). Bioinformatics. 2014 Jun 28;btu410.\n",
	       "\n",
        );
}
