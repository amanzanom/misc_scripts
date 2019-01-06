#!/usr/bin/perl

# Create repeat region annotations from sequence. NOTE: Needs RECOn verwion >=1.07

use Getopt::Long;


GetOptions(\%opts, "infile|i=s", "prefix|p=s", "recon=s", "int:i", "iden:i", "eval:i", "clean!", "help|h!");
        sub usage(){
                die "USAGE :: perl repeatFind.pl -i file -p file -recon RECONPATH [-int integer] [-ident integer] [-eval evalue_string] [-clean] [-help|-h]\n\n
	-infile|i\n\tFasta containing sequences [String]\n
	-prefix|p\n\tPrefix for outfile(s) [String]\n
	-recon\n\tPath to RECON (tested with v1.07) [String]\n
	-int\n\tInteger to use in RECON (see RECON 00README file for details (default 1) [Integer]\n
	-iden\n\tIdentity cut-off for megablast hits (default 97) [Integer]\n
	-eval\n\tE-value cut-off for megablast hits (default 1e-05) [String]\n
	-clean\n\tClean all intermediate files (default 1) [Boolean]\n
	-help|h\n\tGet help for using this script [Boolean]\n";
        }

if ($opts{'help'} || !$opts{'infile'} || !$opts{'prefix'} || !$opts{'recon'}){
	if (!$opts{'infile'}){
		print "infile missing, please check usage manual\n\n";
	}
	if (!$opts{'recon'}){
		print "Root to RECON not specified, please check usage manual\n\n";
	}
	if (!$opts{'prefix'}){
		print "Outfile(s) prefix missing, please check usage manual\n\n";
	}
	&usage;
}

# Variable definition

## Capture options and/or set default values
my $inFile= $opts{'infile'};
my $outDir= $opts{'dir'};
$outDir=~ s/\/+$//;
my $prefix= $opts{'prefix'};
my $reconPath= $opts{'recon'};
$reconPath=~ s/\/+$//;
my $integer=1;
if ($opts{'int'}){
	$integer= $opts{'int'};
}
my $identity=97;
if ($opts{'iden'}){
	$identity= $opts{'iden'};
}
my $eval="1e-05";
if ($opts{'eval'}){
	$eval= $opts{'eval'};
}
my $cleanUp= 0;
if ($opts{'clean'}){
	$cleanUp= $opts{'clean'};
}

## Other variables
my $line='';
my %repeatSeq=();
my @tempArray=();
my $tempElem=0;
my $tempString='';
my $sequence='';
my $i=0;

### Construct seq_name_list_file for RECON
open (GFFTWO, ">$prefix.gff2") || die ("Unable to open file for writing: $prefix.gff2\n$!\n");
open (GFFTHR, ">$prefix.gff3") || die ("Unable to open file for writing: $prefix.gff3\n$!\n");
print GFFTWO "##gff-version 2\n##source-version priScaff v1.0\n";
print GFFTHR "##gff-version 3\n";
print GFFTWO "##source-version RECON v1.07\n";

### Print fasta sequences for GFF v2
$i=0;
open (FASTA, "$inFile") || die ("Unable to open file for reading: $inFile\n$!\n");
while ($line=<FASTA>){
	if ($line=~ m/^>/){
		$line=~ s/^>/DNA /;
		if ($i){
			print GFFTWO "##end-DNA\n";
		}
		$i++;
	}
	$line=~ s/^/##/;
	print GFFTWO $line;
}
close (FASTA);
print GFFTWO "##end-DNA\n";
$i=0;

### Create seq_name_list file for RECON
open (SEQNAMELST, ">$prefix.seqNameLst") || die ("Unable to open file for writting: $prefix.seqNameLst\n$!\n");
print SEQNAMELST `grep -c "^>" $inFile`;
print SEQNAMELST `grep "^>" $inFile | sed 's/^>//' | sort`;
close (SEQNAMELST);

### Run blast, MSPCollect.pl and recon.pl
$cmd="formatdb -i $inFile -p F";
print "$cmd\n";
system ($cmd) == 0 || die "Unable to run formatdb: error code $?\n$!\n";
$cmd="megablast -i $inFile -d $inFile -o $prefix.blast.out -p 99 -W 16 -e 1e-05 -m 8";
print "$cmd\n";
system ($cmd) == 0 || die "Unable to run megablast: error code $?\n$!\n";
$cmd="perl $reconPath/scripts/MSPCollect.pl $prefix.blast.out > $prefix.blast.MSP";
print "$cmd\n";
system ($cmd) == 0 || die "Unable to run MSPCollect.pl: error code $?\n$!\n";
if (-s "$prefix.blast.MSP"){
	$cmd="perl $reconPath/scripts/recon.pl $prefix.seqNameLst $prefix.blast.MSP $integer";
	print "$cmd\n";
	system ($cmd) == 0 || die "Unable to run recon.pl: error code $?\n$!\n";
	### Store repeat positions
	open (ELES, "summary/eles") || die ("Unable to open file for reading: summary/eles\n$!\n");
	while ($line=<ELES>){
		if ($line=~ m/^#/){
			next;
		}
		$line=~ s/^\s+//;
		$line=~ s/\s+$//;
		@tempArray= split(/\s+/, $line);
		$tempString="+";
		if ($tempArray[2] == -1){
			$tempString="-";
		}
		print GFFTWO "$tempArray[3]\tRECON\trepeat_region\t$tempArray[4]\t$tempArray[5]\t.\t$tempString\t.\tlocus_tag=REP_F$tempArray[0]E$tempArray[1]; note=\"family:$tempArray[0] element:$tempArray[1]\"\n";
		print GFFTHR "$tempArray[3]\tRECON\trepeat_region\t$tempArray[4]\t$tempArray[5]\t.\t$tempString\t.\tID=REP_F$tempArray[0]E$tempArray[1];Name=REP_F$tempArray[0]E$tempArray[1];Note=family:$tempArray[0],element:$tempArray[1]\n";
		$repeatSeq{$tempArray[3]}{"REP_F$tempArray[0]E$tempArray[1]"}[0]= $tempArray[4]-1; # -1 to save in 0-based index
		$repeatSeq{$tempArray[3]}{"REP_F$tempArray[0]E$tempArray[1]"}[1]= $tempArray[5]-1; # -1 to save in 0-based index
	}
	close (ELES);
}
$tempArray=();
$tempString='';

close (GFFTWO);

### Print fasta sequences for GFF v3
print GFFTHR "##FASTA\n";
open (FASTA, "$inFile") || die ("Unable to open file for reading: $inFile\n$!\n");
while ($line=<FASTA>){
	print $line;
}
close (FASTA);
print GFFTHR "###\n";
close (GFFTHR);



### Cleanup RECON files
if ($cleanUp){
	$cmd="rm -fr edge_redef_res ele* images summary formatdb.log $inFile.n* $prefix.blast* $prefix.seqNameLst";
	print "$cmd\n";
	system ($cmd) == 0 || die "Unable to run cleanUp of RECON files: error code $?\n$!\n";
}


open (NEWFASTA, ">$prefix.masked.fasta") || die ("Unable to open file for reading: $prefix.masked.fasta\n$!\n")
open (FASTA, "$inFile") || die ("Unable to open file for reading: $inFile\n$!\n");
while ($line=<FASTA>){
	chomp $line;
	if ($line=~ /^>(\S+)/ || eof){
		if (eof){
			$sequence.=uc($line);
		}
		if (length($sequence)){
			print "Masking repeats for scaffold: $ctgId (length=" . length($sequence) . ")\n";
			foreach $tempElem (keys %{$repeatSeq{$ctgId}}){
				$tempString= substr($sequence, $repeatSeq{$ctgId}{$tempElem}[0], ($repeatSeq{$ctgId}{$tempElem}[1]-$repeatSeq{$ctgId}{$tempElem}[0]+1));
				$tempString= lc($tempString);
				print "repeat $tempElem: $repeatSeq{$ctgId}{$tempElem}[0], $repeatSeq{$ctgId}{$tempElem}[1]\n";
				substr($sequence, $repeatSeq{$ctgId}{$tempElem}[0], ($repeatSeq{$ctgId}{$tempElem}[1]-$repeatSeq{$ctgId}{$tempElem}[0]+1), $tempString)
			}
			print NEWFASTA ">$ctgId\n";
			for ($i=0; $i<length($sequence); $i+=70){
				if (length(substr($sequence, $i, length($sequence)-$i))<70){
					print NEWFASTA substr($sequence, $i, length($sequence)-$i)."\n";
				}
				else {
					print NEWFASTA substr($sequence, $i, 70)."\n";
				}
			}
		}
		$ctgId=$1;
		$sequence='';
		next;
	}
	$sequence.=uc($line);
}
close (FASTA);
