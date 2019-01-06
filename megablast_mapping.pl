#!/usr/local/bin/perl

##########################################################
# Perl script to clean BLAST tabbed output by an e-value #
# threshold and/or identity value threshold and filter	 #
# overlapping hits by certain percentage overlap between #
# them (default (50%)					 #
##########################################################


use Getopt::Long;
use Data::Dumper;

GetOptions(\%opts, "infile|i=s", "evalue:f", "identity:f", "overlap:f", "seq:s", "p:s", "rm!", "sum!", "help|h!");
        sub usage(){
                die "USAGE :: cut_out_BLAST_results.pl -infile|i blastout [-evalue FLOAT] [-identity INT] [-overlap FLOAT -seq DBseq.fasta][-rm] [-help|-h]\n
	-infile or -i\n\tFile containing blast results [String]\n
	-evalue\n\tIf present, clean by evalue threshold (minimum) [Float < 10]\n
	-identity\n\tIf present, clean by identity threshold (minimum) [Float < 1]\n
	-overlap\n\tIf present, clean by minimum overlap threshold (minimum) [Integer >=1]\n
	-p\n\tIf proteins compared give name of protein multifasta [String]\n
	-rm\n\tIf present, remove from sequence file\n
	-sum\n\tIf present create a summmary of clean\n
	-help or -h\n\tGet help for using this script\n";
        }

if ($opts{'help'} || !$opts{'infile'} || (!$opts{'evalue'} && !$opts{'identity'} && !$opts{'overlap'})){
	if (!$opts{'infile'} || !$opts{'infile'}){
		print "One or both infiles missing, check usage manual\n";
	}
	if (!$opts{'evalue'} && !$opts{'identity'}){
		print "Both cut thresholds missing, check usage manual\n";
	}
	die &usage();
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

sub getSign {
	my $x=$_[0];
	if ($x==0){
		return 0;
	}
	if ($x<0){
		return -1;
	}
	if ($x>0){
		return 1;
	}
}


#####Main program#####
{
	open (SUM, ">$opts{'infile'}.cut.summary") || die "Unable to open: $opts{'infile'}.cut.summary\n$!";
	if ($opts{'evalue'}){
		print SUM "#E-value threshold: $opts{'evalue'}\n"
	}
	if ($opts{'identity'}){
		print SUM "#Identity threshold: $opts{'identity'}\n"
	}
	if ($opts{'overlap'}){
		print SUM "#Overlap threshold: $opts{'overlap'}\n"
	}
	#####Load sequences from fasta file#####
	if ($opts{'seq'}){
		$flag=0;
		open (IN, "$opts{'seq'}") || die "Unable to open: $opts{'seq'}\n$!";
			while ($line=<IN>){
				chomp $line;
				if ($line=~/^>(\S+)/){
					$name=$1;
					$flag=1;
					$seq{$name}=0;
					next;
				}
				if ($flag){
					$seq{$name}+=length($line);
				}
			}
		close (IN);
		print SUM "#Sequences loaded: ".scalar(keys %seq)."\n";
	}
	
	#####Load blast results#####
	open (IN, "$opts{'infile'}") || die "Unable to open: $opts{'infile'}\n$!";
	open (OUT, ">$opts{'infile'}.cut") || die "Unable to open: $opts{'infile'}.cut\n$!";
		$i=0;
		$c=0;
		print SUM "#Cutoffs\n";
		while ($line=<IN>){
			chomp $line;
			@temp= split (/\s+/, $line);
			$temp[10]=~ s/e/*10**/;
			if ( ($opts{'evalue'} && eval($temp3)>$opts{'evalue'}) || ($opts{'identity'} && $temp[2]<$opts{'identity'}) || ($opts{'overlap'} && $opts{'seq'} && ($temp[3]/$seq{$temp[0]})*100 < $opts{'overlap'}) ){
#			if ( ($opts{'evalue'} && eval($temp4)>$opts{'evalue'}) || ($opts{'identity'} && $temp[2]<$opts{'identity'}) || ($opts{'overlap'} && $opts{'seq'} && $temp[3]<$opts{'overlap'}) ){
				print SUM join("	", @temp)."\n";
				$c++;
			}
			else {
				print OUT "$line\n"
			}
			undef @temp;
			undef $temp2;
			undef $temp3;
		}
	close (IN);
	close (OUT);
	close (SUM);
}

