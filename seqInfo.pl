#!/usr/bin/perl

#===============================================================================
#   Author: Alejandro Manzano Marin
#
#   File: seqInfo.pl
#   Date: 15-11-2012
#   Version: 1.0
#
#   Usage:
#      perl seqInfo.pl -infile|i inFile.fasta [-reads flx|plus] [-prefix|p outFile] [-graph] [options]
#
#      Check out 'perl seqInfo.pl -h' for short usage manual and info on the software.
#
#    Description: This program will analyze a set of FASTA formatted DNA sequences and can give various
#                 statistics as n50, n95, mean G+C content, G+C content of first position, etc.. It can
#                 also make graphs in R to visualize some of the results.
#
#    Contact: Contact the author at alejandro.manzano@uv.es using 'seqInfo: ' as
#             as begining for subject for bug reporting, feature request or whatever
#             (of course realted to the software).
#
#    COPYRIGHT: Copyright (C) 2012  Alejandro Manzano-Marin.
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
use Data::Dumper;
use Pod::Usage;
use Cwd 'abs_path';


# Define subroutines
sub version {
	my $software= $_[0];
	my $version= $_[1];
	print "$software v$version\n";
	exit (0);
}

sub countGC {
	my $seq=$_[0];
	my $gcCount=0;
	while ($seq=~ m/[GC]/gi){
		$gcCount++;
	}
	return ($gcCount);
}

sub extractPos {
	my $seq=$_[0];
	my $pos1=$_[1];
	my $pos2=$_[2];
	my $pos3=$_[3];
	my $newSeq='';
	my $tempCodon='';
	while ($seq=~ m/(\w{1,3})/gi){
		$tempCodon=$1;
		if ($pos1){
			$newSeq.=substr($tempCodon, 0, 1);
		}
		if ($pos2){
			$newSeq.=substr($tempCodon, 1, 1);
		}
		if ($pos3){
			$newSeq.=substr($tempCodon, 2, 1);
		}
	}
	return ($newSeq);
}

sub sum {
	my $refArray=$_[0];
	my $i=0;
	my $total=0;
	for ($i=0; $i<scalar(@{$refArray}); $i++){
		$total+=$refArray->[$i];
	}
	return ($total);
}

sub stdev {
	my $refArray=$_[0];
	my $mean=&sum($refArray)/scalar(@{$refArray});
	my $i=0;
	my $total=0;
	for ($i=0; $i<scalar(@{$refArray}); $i++){
		$total+=(($refArray->[$i])-$mean)**2;
	}
	$total/=scalar(@{$refArray});
	return (sqrt($total));
}

sub max {
	my $refArray=$_[0];
	my $max=$refArray->[0];
	my $i=1;
	for ($i=1; $i<scalar(@{$refArray}); $i++){
		if ($refArray->[$i]>$max){
			$max=$refArray->[$i];
		}
	}
	return ($max);
}

sub min {
	my $refArray=$_[0];
	my $min=$refArray->[0];
	my $i=1;
	for ($i=1; $i<scalar(@{$refArray}); $i++){
		if ($refArray->[$i]<$min){
			$min=$refArray->[$i];
		}
	}
	return ($min);
}

sub Nx {
	my $refArray=$_[0];
	my $total=$_[1];
	my $x=$_[2];
	my $i=0;
	my $sum=0;
	my $Nx=0;
	my @tempArray=sort{ $b <=> $a } @{$refArray};
	my $x= $x/100;
	for (my $i=0;$i<@tempArray;$i++){
		$sum+=$tempArray[$i];
		if ($sum>=$total*$x){
			$Nx=$tempArray[$i];
			last;
		}
	}
	undef(@tempArray);
	return($Nx);
}

sub Lx {
	my $refArray=$_[0];
	my $total=$_[1];
	my $x=$_[2];
	my $i=0;
	my $sum=0;
	my $Lx=0;
	my @tempArray=sort{ $b <=> $a } @{$refArray};
	my $x= $x/100;
	for (my $i=0;$i<@tempArray;$i++){
		$sum+=$tempArray[$i];
		if ($sum>=$total*$x){
			$Lx=$i+1;
			last;
		}
	}
	undef(@tempArray);
	return($Lx);
}

# Variable definition

## Define other variables
my $scriptPath=abs_path($0);
$scriptPath=~ s/[^\/]+$//;

my @seqLengths=();
my $line='';
my $length=0;
my $seq='';
my $id='';
my $contig=0;
my $i= 0; # counter
my %temp= ('num' => 0,
	'seqLength' => 0);
my $total_length= 0;

## General variables
my $PROGRAMNAME= "seqInfo";
my $VERSION= "1.0";

## Define options default values
my $inFile= "";
my $reads='0';
my $ignore_n=0;

my $prefix= "seqInfo_out";
my $graph= 0;
my $verbose= 0;

my @Nx= ();
my @Lx= ();
my $calcGC=0;
my $calc1GC=0;
my $calc2GC=0;
my $calc12GC=0;
my $calc3GC=0;
my $all=0;

## Define options hash
GetOptions(\%opts, 
	"infile|i=s", 
	"reads:s", 
	"ignoren!", 
	"prefix|p:s", 
	"graph!", 
	"Nx:i@", 
	"Lx:i@", 
	"GC!", 
	"1GC!", 
	"2GC!", 
	"12GC!", 
	"3GC!", 
	"all!", 
	"man!", 
	"help|h!") || pod2usage(2);

if ($opts{'help'}){
	pod2usage(-exitval => 1,  -verbose => 1);
}
if ($opts{'man'}){
	pod2usage(-exitval => 0, -verbose => 2);
}

## Script documetation
=pod

=head1 NAME

seqInfo

=head1 VERSION

seqInfo v1.0

=head1 SYNOPSIS

perl seqInfo.pl --infile|-i inFile.fasta [--reads flx|plus] [--prefix|-p outFile] [--graph] [--Nx integer [--Nx integer ...]] [--Lx integer [--Lx integer ...]] [--GC] [--1GC] [--2GC] [--12GC] [--3GC] [--all] [--verbose|-v] [--help|-h] [--man] [--version]

=head1 DESCRIPTION

This program will analyze a set of FASTA formatted DNA sequences and can give various statistics as n50, n95, mean G+C content, G+C content of first position, etc.. It can also make graphs in R to visualize some of the results.

=head1 OPTIONS

=head2 INPUT

=over 8

=item B<-i> | B<--infile> <string> (mandatory)

File containing sequences to be analyzed in FASTA format.

=item B<--reads> <string> (default: 0)

Indicates that input reads are sequencing reads (currently only 454) of either "flx" or "plus" and graphs information useful for this type of data.

=item B<--ignoren> <string> (default: 0)

If used, the script will ignore undefined nucleotides ('N' or 'n', normally arising from assembly gaps) for the calculation of the stats.

=back

=head2 OUTPUT

=over 8

=item B<-p> | B<--prefix> <string> (default: "seqInfo_out")

Prefix used for outFiles.

=item B<--graph> <boolean> (default: 0)

If used, make graphs using R.

=item B<-v> | B<--verbose> <boolean> (default: 0)

Prints status and info messages during processing.

=back

=head2 ANALYSES REPORT OPTIONS

=over 8

=item B<--Nx> <integer>

Calculate Nx stats (makes more sense for example with assembly contigs or scaffolds). Can be used more than once for many Nx values (e.g. N50, N80).

=item B<--Lx> <integer>

Calculate Lx stats (makes more sense for example with assembly contigs or scaffolds). Can be used more than once for many Lx values (e.g. L50, L80).

=item B<--GC> <boolean>

alculate G+C mean content of sequence(s).

=item B<--1GC> <boolean>

Calculate G+C of first positions of sequence(s) (Makes sense if analyzing CDSs).

=item B<--2GC> <boolean>

Calculate G+C of second positions of sequence(s) (Makes sense if analyzing CDSs).

=item B<--12GC> <boolean>

Calculate G+C of first+second positions of sequence(s) (Makes sense if analyzing CDSs).

=item B<--3GC> <boolean>

Calculate G+C of third positions of sequence(s) (Makes sense if analyzing CDSs).

=item B<--all> <boolean>

Calculate all stats.

=back

=head2 INFO AND HELP

=over 8

=item B<-h> | B<--help> <boolean>

Print useful help on using this script.

=item B<--man> <boolean>

Print the full documentation.

=item B<--version> <boolean>

Print program version.

=back

=head1 AUTHOR

Alejandro Manzano-Marin, C<< <alejandro_dot_manzano_at_uv_dot_es> >>

=head1 BUGS

If you find a bug please email me at C<< <alejandro_dot_manzano_at_uv_dot_es> >> so that I can keep making seqInfo better.

=head1 COPYRIGHT

Copyright (C) 2012  Alejandro Manzano-Marin.

=head1 LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut

# Check necessary options and paths
if ($opts{'help'} || !$opts{'infile'} || ($opts{'reads'} && $opts{'reads'}!~ m/^(flx|plus)$/)){
	if (!$opts{'help'}){
		print STDERR "ERROR:\n";
		if (!$opts{'infile'}){
			print STDERR "FASTA infile missing\n";
		}
		if ($opts{'reads'}!~ m/^(flx|plus)$/){
			print STDERR "Reads sequencing platform not recognized, please use 'flx' or 'plus'\n";
		}
		print STDERR "Please check usage manual\n";
	}
	print STDERR "\n";
	pod2usage(-exitval => 1,  -verbose => 1);
}

## Capture options
$inFile= $opts{'infile'};
if ($opts{'reads'}){
	$reads= $opts{'reads'};
}
if ($opts{'ignoren'}){
	$ignore_n= $opts{'ignoren'};
}

if ($opts{'prefix'}){
	$prefix= $opts{'prefix'};
}
if ($opts{'graph'}){
	$graph=1;
}
if ($opts{'verbose'}){
	$verbose= 1;
}


if ($opts{'Nx'} || $all){
	if ($all){
		$Nx[0]=50;
	}
	if ($opts{'Nx'}){
		push(@Nx, @{$opts{'Nx'}});
	}
}
if ($opts{'Lx'} || $all){
	if ($all){
		$Lx[0]=50;
	}
	if ($opts{'Lx'}){
		push(@Lx, @{$opts{'Lx'}});
	}
}
if ($opts{'GC'} || $all){
	$calcGC=1;
	my @seqGC=();
}
if ($opts{'1GC'} || $all){
	$calc1GC=1;
	my @seqGC1st=();
}
if ($opts{'2GC'} || $all){
	$calc2GC=1;
	my @seqGC2nd=();
}
if ($opts{'12GC'} || $all){
	$calc12GC=1;
	my @seqGC1st2nd=();
}
if ($opts{'3GC'} || $all){
	$calc3GC=1;
	my @seqGC3rd=();
}
if ($opts{'all'}){
	$all=1;
}



if ($opts{'version'}){
	&version($PROGRAMNAME, $VERSION);
}


# Main program
open (AUX , ">$prefix.aux") or die ("Unable to create $prefix.aux\n$!\n");
open (FASTA, "<$inFile") or die ("Unable to open $inFile\n$!\n");
while ($line=<FASTA>){
	chomp ($line);
	if ($line =~ /^>(\S+)/ || eof){
		if (eof){
			$seq .= $line;
		}
		if (length($seq)>0){
			if ($ignore_n){
				$seq=~ s/n+//gi;
			}
			$temp{'seqLength'}= length($seq);
			push(@seqLengths, $temp{'seqLength'});
			print AUX $id . "\t" . $temp{'seqLength'} . "\t";
			if ($calcGC || $all){
				$temp{'num'}= &countGC($seq);
				push(@seqGC, $temp{'num'});
				print AUX ($temp{'num'}/$temp{'seqLength'}) . "\t";
			}
			else {
				print AUX "NULL\t";
			}
			if ($calc1GC || $all){
				$temp{'num'}= &countGC(&extractPos($seq, 1, 0, 0));
				$temp{'seqLength'}= length(&extractPos($seq, 1, 0, 0));
				push(@seqGC1st, $temp{'num'}/$temp{'seqLength'});
				print AUX ($temp{'num'}/$temp{'seqLength'}) . "\t";
			}
			else {
				print AUX "NULL\t";
			}
			if ($calc2GC || $all){
				$temp{'num'}= &countGC(&extractPos($seq, 0, 1, 0));
				$temp{'seqLength'}= length(&extractPos($seq, 0, 1, 0));
				push(@seqGC2nd, $temp{'num'}/$temp{'seqLength'});
				print AUX ($temp{'num'}/$temp{'seqLength'}) . "\t";
			}
			else {
				print AUX "NULL\t";
			}
			if ($calc12GC || $all){
				$temp{'num'}= &countGC(&extractPos($seq, 1, 1, 0));
				$temp{'seqLength'}= length(&extractPos($seq, 1, 1, 0));
				push(@seqGC1st2nd, $temp{'num'}/$temp{'seqLength'});
				print AUX ($temp{'num'}/$temp{'seqLength'}) . "\t";
			}
			else {
				print AUX "NULL\t";
			}
			if ($calc3GC || $all){
				$temp{'num'}= &countGC(&extractPos($seq, 0, 0, 1));
				$temp{'seqLength'}= length(&extractPos($seq, 0, 0, 1));
				push(@seqGC3rd, $temp{'num'}/$temp{'seqLength'});
				print AUX ($temp{'num'}/$temp{'seqLength'}) . "\n";
			}
			else {
				print AUX "NULL\n";
			}
		}
		if (eof){
			last;
		}
		$id=$1;
		$contig++;
		$seq= '';
	}
	else {
		$seq = $seq . $line;
	}
}
close (FASTA);
close (AUX);



## Print sequences summary

$total_length= &sum(\@seqLengths);
open (SUM , ">$prefix.sum") or die ("Unable to create $prefix.sum\n$!\n");
print SUM "# Sequence Info\n\n";
print SUM "# Total sequences captured:\t$contig\n";
print SUM "# Total bases captured:\t" . $total_length . "\n";
if (length(@Nx) > 0 || $all){
	for ($i= 0; $i<scalar(@Nx); $i++) {
		print SUM "# N". $Nx[$i] . " of sequences captured:\t". &Nx(\@seqLengths, $total_length, $Nx[$i]) ."\n";
	}
}
if (length(@Lx) > 0 || $all){
	for ($i= 0; $i<scalar(@Lx); $i++) {
		print SUM "# L". $Lx[$i] . " of sequences captured:\t". &Lx(\@seqLengths, $total_length, $Lx[$i]) ."\n";
	}
}
print SUM "# Length of smallest sequence (bp):\t" . &min(\@seqLengths) . "\n";
print SUM "# Length of biggest sequence (bp):\t" . &max(\@seqLengths) . "\n";
print SUM "# Average Length of sequences(bp):\t". sprintf("%.4f", ($total_length/scalar(@seqLengths))) ."\n";
print SUM "# StdDev Length of sequences (bp):\t". sprintf("%.4f", &stdev(\@seqLengths)) ."\n";
if ($calcGC || $all){
	print SUM "# Total G+C content:\t". sprintf("%.4f", (&sum(\@seqGC)/$total_length)) ."\n";
	$temp{'num'}= 0;
	for ($i= 0; $i<scalar(@seqGC); $i++) {
		$temp{'num'}+= ($seqGC[$i]/$seqLengths[$i]);
	}
	print SUM "# Mean G+C content:\t". sprintf("%.4f", ($temp{'num'}/scalar(@seqLengths))) ."\n";
}
if ($calc1GC || $all){
	print SUM "# 1st pos. mean G+C content:\t" . sprintf("%.4f", (&sum(\@seqGC1st)/scalar(@seqGC1st))) . "\n";
}
if ($calc2GC || $all){
	print SUM "# 2nd pos. mean G+C content:\t" . sprintf("%.4f", (&sum(\@seqGC2nd)/scalar(@seqGC2nd))) . "\n";
}
if ($calc12GC || $all){
	print SUM "# 1st+2nd pos. mean G+C content:\t" . sprintf("%.4f", (&sum(\@seqGC1st2nd)/scalar(@seqGC1st2nd))) . "\n";
}
if ($calc3GC || $all){
	print SUM "# 3rd pos. mean G+C content:\t" . sprintf("%.4f", (&sum(\@seqGC3rd)/scalar(@seqGC3rd))) . "\n";
}
close (SUM);

## Write R script and execute to graph if appointed by the user by the -graph option
if ($graph){
	$cmd="R --slave --vanilla --args " . $prefix . ".aux " . $prefix . " " . $reads . " < " . $scriptPath . "seqInfo.R";
	system ($cmd);
}


exit (0);
