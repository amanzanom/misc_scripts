#!/usr/bin/perl

use Getopt::Long;
use Data::Dumper;
use Cwd 'abs_path';


GetOptions(\%opts, "infile|i=s", "backbone|b=s", "minOvlp=s", "keepEnds!", "help|h!");
        sub usage(){
                die "USAGE :: perl cleanSamBWT.pl -i inFile -b backboneFile.fasta -minOvlp percentage100Scale [-help|-h]\n\n
	-infile|i\n\tSam file of reads vs the backbone [String]\n\n
	-backbone|b\n\tFASTA of backbone sequence the backbone [String]\n\n
	-minOvlp\n\tMinimum alignment overlap percentage to keep [Real, (0,100] ]\n\n
	-keepEnds\n\tKeep alignments shorter than minAln length if they are at the ends (make sense for expanding backbones) [Boolean]\n\n
	-help or -help\n\tGet help for using this script\n\n";
        }

if ($opts{'help'} || !$opts{'infile'} || !$opts{'minOvlp'}){
	if (!$opts{'infile'}){
		print "infile missing, please check usage manual\n\n";
	}
	if (!$opts{'minOvlp'}){
		print "No minimum alignement length specified, won't run the script cause it doesn't make sense, sorry, please check usage manual\n\n";
	}
	&usage;
}

# Variable definition

## Capture options
my $inFile= $opts{'infile'};
my $backboneFile= $opts{'backbone'};
my $minOverlap= $opts{'minOvlp'};
my $keepEnds= $opts{'keepEnds'};

## Other variables
my %backboneLength=();
my $line='';
my $cigarString='';
my $start=0;
my $strand=0;
my $sequence='';
my $id='';
my $skippedLeft=0;
my $skippedRight=0;


# Main program

open (BACKBONE, "$backboneFile") or die ("Unable to open file for reading: $backboneFile\n$!\n");
while ($line=<BACKBONE>){
	chomp ($line);
	if ($line =~ /^>(\S+)/){
		$id=$1;
		$backboneLength{$id}= 0;
		next;
	}
	$backboneLength{$id}+= length($line);
}
close (BACKBONE);

open (OUTSAM , ">$inFile.clean") or die ("Unable to create file for writting: $inFile.clean\n$!\n");
open (INSAM, "$inFile") or die ("Unable to open file for reading: $inFile\n$!\n");
while ($line=<INSAM>){
	if ($line =~ /^@/){
		next;
	}
	if ($line=~ m/^\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+\S+\s+\S+\s+(\S+)/){
		$strand=$1;
		$id=$2;
		$start=$3;
		$cigarString=$4;
		$sequence=$5;
		$alnLength=0;
		while ($cigarString=~ m/(\d+)[MI]/g){
			$alnLength+=$1;
		}
		$skippedLeft=0;
		$skippedRight=0;
		if ($cigarString=~ m/^(\d+)S/){
			$skippedLeft=$1;
		}
		if ($cigarString=~ m/(\d+)S$/){
			$skippedRight=$1;
		}
		if ($keepEnds && ($alnLength/length($sequence))*100 >= $minOverlap || (($alnLength/length($sequence))*100 < $minOverlap && ( ($start - $skippedLeft <= 0 && $skippedRight<=4)  || ($start+$skippedLeft+$alnLength+$skippedRight >= $backboneLength{$id} && $skippedLeft<=4) ) )){
			print OUTSAM $line;
		}
	}
}
close (INSAM);
close (OUTSAM);
