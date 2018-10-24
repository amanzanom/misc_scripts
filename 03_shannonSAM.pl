#!/usr/bin/perl

#===============================================================================
#   Author: Alejandro Manzano-Marin,
#
#   File: shannonSAM
#   Date: 01-04-2013
#   Version: 0.1
#
#   Usage:
#      perl shannonSAM.pl -infile|i infile.sam [-infile|i infile.sam ...] -backbone|b backbone.fasta [-seqType 'nuc'|'aa'] [options]
#
#      Check out 'perl shannonSAM.pl -h' for short usage manual and info on the software.
#
#    Description: shannonSAM is a script designed to calculate the relative frequency of A,T,G,C and - characters
#                 to calculate Shannon entropy for positions in a SAM alignment considering A,T,G,C and - using
#                 the script plotShannonSeq.R. The calculations are done as described in Zagordi et al Read length
#                 vrsus Depth Coverage for Viral Quasispecies Recnstruction. PLoS ONE 7(10) e47046.
#
#    Contact: Contact the author at alejandro.manzano@uv.es using 'shannonSAM: ' as
#             as begining for subject for bug reporting, feature request or whatever
#             (of course realted to the software).
#
#    COPYRIGHT: Copyright (C) 2013  Alejandro Manzano-Marin.
#
#    LICENCE: This program is free software: you can redistribute it and/or modify it under the terms
#             of the GNU General Public License as published by the Free Software Foundation, either
#             version 3 of the License, or (at your option) any later version.
#             This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
#             without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#             See the GNU General Public License for more details.
#             You should have received a copy of the GNU General Public License along with this program.
#             If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================


# Load modules
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;


# Define subroutines
sub printVersion {
	my $software= $_[0];
	my $version= $_[1];
	my $fileHandle= $_[2];
	print $fileHandle "$software v$version\n";
	exit (0);
}

sub max {
	my @array=@_;
	my $max=$array[0];
	for (my $i=1; $i<scalar(@array); $i++){
		if ($array[$i]>$max){
			$max=$array[$i];
		}
	}
	return ($max);
}

sub min {
	my @array=@_;
	my $min=$array[0];
	for (my $i=1; $i<scalar(@array); $i++){
		if ($array[$i]<$min){
			$min=$array[$i];
		}
	}
	return ($min);
}

sub seqOvlpSize {
	my ($x1, $x2, $y1, $y2)=@_;
	my $overlap=0;
	if ((&min($x1,$x2)-&max($y1,$y2))*(&max($x1,$x2)-&min($y1,$y2)) > 0){
		return (0);
	}
	my @temp= sort { my $a <=> my $b } ($x1, $x2, $y1, $y2);
#	@temp= sort { $a <=> $b } ($x1, $x2, $y1, $y2);
	pop @temp;
	shift @temp;
	$overlap= abs($temp[0]-$temp[1])+1;
	return ($overlap);
}

# Variable definition

## Define other variables
my @abc= ();
my $line= '';
my $id= '';
my %backboneSeq= ();
my $readID= '';
my $refID= '';
my $alnStart= 0;
my $alnCigar= '';
my $readSeq= '';
my $alnLength= 0;
my %nucAbsFreq= ();
my %temp= ('var' => '', 
	'num' => 0, 
	'seq' => '', 
	'length' => 0, 
	'start' => 0);
my $i= 0;
my $file= '';
my %flag= ('print' => 0);

## General variables
my $PROGRAMNAME= 'shannonSAM';
my $VERSION= '0.1';


## Define options default values
my @opt_inFiles= ();
my $opt_backboneFile= '';

my $opt_seqType= 'nuc';

my $opt_verbose= 0;
my $opt_man= 0;
my $opt_help= 0;
my $opt_printVersion= 0;


## Define options hash
GetOptions(\%opts, 
	'infile|i=s@' => \@opt_inFiles, 
	'backbone|b=s' => \$opt_backboneFile, 
	'seqType:s' => \$opt_seqType, 
	'verbose|v!' => \$opt_verbose, 
	'help|h!' => \$opt_help, 
	'man!'  => \$opt_man, 
	'version!' => \$opt_printVersion) || pod2usage(-exitval => 1,  -verbose => 2);
	
if ($opt_help){
	pod2usage(-exitval => 1,  -verbose => 1);
}
if ($opt_man){
	pod2usage(-exitval => 0, -verbose => 2);
}
if ($opt_printVersion){
	&printVersion($PROGRAMNAME, $VERSION, \*STDERR);
}
	

# Script documetation

=pod

=head1 NAME

shannonSAM - Calculate Shannon entropy for an alignent in SAM format.

=head1 VERSION

shannonSAM v0.1

=head1 SYNOPSIS

perl shannonSAM.pl -infile|i infile.sam [-infile|i infile.sam ...] -backbone|b backbone.fasta [-seqType 'nuc'|'aa'] [-range "start-end"] [-v|-verbose] [-h|-help] [-man] [-version]

=head1 DESCRIPTION

shannonSAM is a script designed to extract reads from a certain window of a SAM alignment file.

=head1 OPTIONS

=head2 INPUT

=over 8

=item B<-infile> | B<-i> <string> (mandatory)

SAM alignment file(s) (tested using bowtie2 v2.1.0 SAM output).

=item B<-backbone> | B<-b> <string> (mandatory)

Backbone used for mapping in FASTA format.

=back

=head2 SETTINGS

=over 8

=item B<-seqType> <string> (default: 'nuc')

Type of sequence, either 'nuc' for nucleotides or 'aa' for aminoacids.

=back

=head2 INFO AND HELP

=over 8

=item B<-v> | B<-verbose> <boolean>

Prints status and info messages while processing.

=item B<-h> | B<-help> <boolean>

Print useful help on using this script.

=item B<-man> <boolean>

Print the full documentation.

=item B<-version> <boolean>

Print program version.

=back

=head1 AUTHOR

Alejandro Manzano-Marin, C<< <alejandro_dot_manzano_at_uv_dot_es> >>

=head1 BUGS

If you find a bug please email me at C<< <alejandro_dot_manzano_at_uv_dot_es> >> so that I can keep making fastx_filter better.

=head1 COPYRIGHT

Copyright (C) 2012  Alejandro Manzano-Marin.

=head1 LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut


# Assign values to variables dependant on options


# Check necessary options and paths
if (!@opt_inFiles || !$opt_backboneFile || ($opt_seqType && $opt_seqType!~ m/^(nuc|aa)$/)){
	print STDERR "ERROR:\n";
	if (!@opt_inFiles){
		print STDERR "SAM alignment input file(s) missing\n";
	}
	if (!$opt_backboneFile){
		print STDERR "No backbone file specified\n";
	}
	if ($opt_seqType && $opt_seqType!~ m/^(nuc|aa)$/){
		print STDERR "Range is not correctly specified\n";
	}
	print STDERR "Please check usage manual\n";
	pod2usage(-exitval => 1,  -verbose => 1);
}


# Main program

## Define abc
#if ($opt_seqType eq 'nuc'){
#	push(@abc,'A','T','G','C','GAP');
#}
#elsif ($opt_seqType eq 'aa'){
##	push();
#}
%abc= (
	'nuc' => {
		'A' => 'A',
		'T' => 'T',
		'G' => 'G',
		'C' => 'C',
		'-' => 'GAP'
	},
	'aa' => {
		'A' => 'A',
		'B' => 'B',
		'-' => 'GAP'
	}
);
	
#print Dumper %abc;$pene= <STDIN>;

## Read backbone file 
if ($opt_verbose){
	print STDERR "Reading backbone reference file $backboneFile.. ";
}

open (BCKBIN, "<$opt_backboneFile") || die ("Unable to open file for reading: $opt_backboneFile\n$!\n");
while ($line=<BCKBIN>){
	chomp $line;
	if ($line=~ m/^>(\S+)/){
		$refID= $1;
		$backboneSeq{$refID}= '';
		next;
	}
	$backboneSeq{$refID}.= $line;
}
close (BCKBIN);

if ($opt_verbose){
	print STDERR "DONE\n";
}


#foreach $refID (keys %backboneSeq){
#	@{$nucAbsFreq{$refID}}= ();
#	for ($i=0; $i<length($backboneSeq{$refID}); $i++){
#		%{$nucAbsFreq{$refID}[$i]}= ( 'A' => 0, 'T' => 0, 'G' => 0, 'C' => 0, '-' => 0);
#	}
#}

## Read SAM file and store A,T,G,C or - characters
foreach $file (@opt_inFiles){
	if ($opt_verbose){
		print STDERR "Reading SAM file $file.. ";
	}
	
	open (SAMIN, "<$file") || die "Unable to open file $file for reading\n$!\n";
	while ($line= <SAMIN>){
#		print $line;$pene=<STDIN>;
		chomp $line;
		if ($line=~ m/^@/){
			next;
		}
		if ($line=~ m/^(\S+)\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+\S+\s+\S+\s+(\S+)/){
			$readID= $1;
			$refID= $2;
			$alnStart= $3;
			$alnCigar= $4;
			$readSeq= $5;
			$alnLength= 0;
			if (!defined($nucAbsFreq{$refID})){
				@{$nucAbsFreq{$refID}}= ();
			}
			$temp{'start'}= 0;
			$temp{'seq'}= '';
			while ($alnCigar=~ m/(\d+[MIDSH])/g){ # Go through the cigar string to remove insertions
				$temp{'var'}= $1;
				if ($temp{'var'}=~ m/(\d+)M/){
					$temp{'length'}= $1;
					$temp{'seq'}.= substr($readSeq, $temp{'start'}, $temp{"length"});
					$temp{'start'}+= $temp{'length'};
				}
				elsif ($temp{'var'}=~ m/(\d+)I/){
					$temp{'length'}= $1;
					$temp{'start'}+= $temp{'length'};
				}
				elsif ($temp{'var'}=~ m/(\d+)D/){
					$temp{'length'}= $1;
					$temp{'seq'}.= '-';
				}
				elsif ($temp{'var'}=~ m/(\d+)[SH]/){
					$temp{'length'}= $1;
					$temp{'start'}+= $temp{'length'};
				}
			}
			if (defined($nucAbsFreq{$refID})){
				for ($i=0; $i<length($temp{'seq'}); $i++){
					$nucAbsFreq{$refID}[$alnStart+$i-1]{$abc{$opt_seqType}{uc(substr($temp{'seq'}, $i, 1))}}++;
				}
			}
		}
	}
	close (SAMIN);
	
	if ($opt_verbose){
		print STDERR "DONE\n";
	}
}

## Print file of character relative frequencies
if ($opt_verbose){
	print STDERR "Printing relative frequency file.. ";
}

foreach $refID (sort {lc($a) cmp lc($b)} keys %nucAbsFreq){
	print STDOUT '## sequence ' . $refID . "\n";
	print STDOUT 'pos';
	foreach $letter (sort {lc($a) cmp lc($b)} keys %{$abc{$opt_seqType}}){
		print STDOUT "\t" . $abc{$opt_seqType}{$letter};
	}
	print STDOUT "\n";
	for ($i=0; $i<length($backboneSeq{$refID}); $i++){
		print STDOUT ($i+1);
		foreach $letter (sort {lc($a) cmp lc($b)} keys %{$abc{$opt_seqType}}){
			if (defined($nucAbsFreq{$refID}[$i]{$abc{$opt_seqType}{$letter}})){
				print STDOUT "\t" . $nucAbsFreq{$refID}[$i]{$abc{$opt_seqType}{$letter}};
			}
			else {
				print STDOUT "\t0";
			}
		}
		print STDOUT "\n";
#		print STDOUT join("	", $nucAbsFreq{$refID}[$i]{'A'}, $nucAbsFreq{$refID}[$i]{'T'}, $nucAbsFreq{$refID}[$i]{'G'}, $nucAbsFreq{$refID}[$i]{'C'}, $nucAbsFreq{$refID}[$i]{'-'}) . "\n";
#		print STDOUT join("	", $nucAbsFreq{$refID}[$i]{'A'}, $nucAbsFreq{$refID}[$i]{'T'}, $nucAbsFreq{$refID}[$i]{'G'}, $nucAbsFreq{$refID}[$i]{'C'}, $nucAbsFreq{$refID}[$i]{'-'}) . "\n";
	}
}

if ($opt_verbose){
	print STDERR "DONE\n";
}


exit (0);
