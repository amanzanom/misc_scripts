#!/usr/bin/perl

use Getopt::Long;


GetOptions(\%opts, "infile|i=s", "prefix|p=s", "seqType=s", "graph!", "help|h!");
        sub usage(){
                die "USAGE :: perl seqInfo.pl -i file -p file [-first] [-second] [-third] [-help|-h]\n\n
	-infile|i\n\tFasta containing sequences [String]\n
	-prefix|p\n\tPrefix for outfile(s) [String]\n
	-seqType\n\tType of DNA sequence, either chromosome, plasmid or viral [String]\n
	-graph\n\tComma separated list of formats to export graph with (pdf and/or ps; be sure to have R installed and in path) [String]\n
	-help|h\n\tGet help for using this script [Boolean]\n";
        }

if ($opts{'help'} || !$opts{'infile'} || !$opts{'prefix'}){
	if (!$opts{'infile'}){
		print "infile missing, please check usage manual\n\n";
	}
	if (!$opts{'prefix'}){
		print "Outfile(s) prefix missing, please check usage manual\n\n";
	}
	&usage;
}

# Check prerequisite software (Defaults were done according to a local install)
#if (!(`which originx`)){
#	print "originx missing or not in PATH\n\n";
#	&usage;
#}
if (!(`which R`)){
	print "R missing or not in PATH\n\n";
	&usage;
}


# Variable definition

## BIN location variables
$ORIXBIN= "/home/manzanoa/software/ANNOTATION/OTHER/originX/bin/originx";

## Capture options and/or set default values
my $inFile= $opts{'infile'};
my $prefix= $opts{'prefix'};
my $seqType= $opts{'seqType'};
my $orixOpt= "";
my $wSize= "";
if ($seqType eq "chromosome"){ # entire chromosomes as described in Peder Worning et al. Origin of replication in circular prokaryotic genomes. Environmental microbiology(2006)
	$orixOpt= "-o 8 -r 5 -s 1000";
	$wSize= "20,50,55,60,65,70";
}
if ($seqType eq "viral" || $seqType eq "plasmid"){ # viral genome or plasmid sequences as described in Peder Worning et al. Origin of replication in circular prokaryotic genomes. Environmental microbiology(2006)
#	$orixOpt= "-o 8 -r 5 -s 100";
	$orixOpt= "-o 8 -r 5 -s 10";
	$wSize= "20,50,55,60,65,70";
}
### originx options explanation
### -m use median percent instead of percent
### -p percentage of total genome size to use as window for search (default 90)
### -ps ??????
### -o oligolength to use fot oligonucleotide skew (default 6)
### -r regulator (????; default 1.000)
### -t threshold (?????; default 40.000)
### -s step of hypothetical origin positions (default 1000)
### infosum= strand bias

my $graph= $opts{'graph'};

## Check files and directories
if (!(-e $inFile)){
	print "Input file doesn't exist or not found please check\n";
	&usage;
}

## Other variables
my $line= '';
my $w = 0;
my $cmd= '';


# Main program

# Run originx
$wSize=~ s/,$//;
foreach $w (split (/\,/, $wSize)){
	$cmd="$ORIXBIN -p $w $orixOpt $inFile > $prefix.$w.orix";
	print "--> Running originxwith window size $w: $cmd ..";
	system ($cmd) == 0 || die "Unable to run originx: error code $?\n$!\n";
	print "DONE\n";
}

# Plot
if ($graph){
	print "--> Creating R script ..";
	open (TMPRSCRIPT, ">$prefix.temp.R") || die ("Unable to open file for writing: $prefix.temp.R\n$!\n");
	foreach $w  (split (/\,/, $wSize)){
		print TMPRSCRIPT "# Read originx output for window size $w\n";
		print TMPRSCRIPT "skewTable$w <- read.table(file=\"$prefix.$w.orix\", header=FALSE, sep=\"\ \", col.names=c(\"Position\", \"Infosum\", \"GC\", \"AT\"), colClasses=\"numeric\", comment.char = \"#\");\n";
		print TMPRSCRIPT "infoSum$w <- sum(skewTable$w\$Infosum)\n";
		print TMPRSCRIPT "gcWMsum$w <- sum(skewTable$w\$GC)\n";
		print TMPRSCRIPT "atWMsum$w <- sum(skewTable$w\$AT)\n\n";
	}
	print TMPRSCRIPT "# Make vectors to graph\n";
	print TMPRSCRIPT "position <- (skewTable20\$Position)/1000\n";
	print TMPRSCRIPT "strandDiffMedian <- apply(cbind(skewTable50\$Infosum*(infoSum60/infoSum50), skewTable55\$Infosum*(infoSum60/infoSum55), skewTable60\$Infosum, skewTable65\$Infosum*(infoSum60/infoSum65), skewTable70\$Infosum*(infoSum60/infoSum70)), 1 , median)/1000\n";
	print TMPRSCRIPT "strandDiff20 <- (skewTable20\$Infosum)/1000\n";
	print TMPRSCRIPT "gcWeightedMedian <- apply(cbind(skewTable50\$GC, skewTable55\$GC, skewTable60\$GC, skewTable65\$GC, skewTable70\$GC), 1 , median)/1000\n";
	print TMPRSCRIPT "atWeightedMedian <- apply(cbind(skewTable50\$AT, skewTable55\$AT, skewTable60\$AT, skewTable65\$AT, skewTable70\$AT), 1 , median)/1000\n";
	print TMPRSCRIPT "\n# Start graph pdf output\n";
	print TMPRSCRIPT "pdf(\"$prefix.pdf\")\n";
	print TMPRSCRIPT "\tlegend<- c(\"Strand bias, median\", \"Strand bias, 20%\", \"G/C weigthed median\", \"A/T weighted median\")\n";
	print TMPRSCRIPT "\tlegendColor<- c(\"red\", \"darkorange\", \"forestgreen\",\"royalblue4\")\n";
	print TMPRSCRIPT "\tyMax <- max(abs(c(strandDiffMedian, strandDiff20, gcWeightedMedian, atWeightedMedian)))\n";
	print TMPRSCRIPT "\tplot(1, type=\"n\", axes=TRUE, xlim=c(0, max(position)), ylim=c(yMax*-1, yMax), bty=\"o\", xlab=\"position in the gneome\", ylab=\"strand bias/1000\", main=\"$prefix\");\n";
	print TMPRSCRIPT "\tlines(x=position, y=strandDiffMedian, col=\"red\", lwd=2);\n";
	print TMPRSCRIPT "\tlines(x=position, y=strandDiff20, col=\"darkorange\", lwd=2);\n";
	print TMPRSCRIPT "\tlines(x=position, y=gcWeightedMedian, col=\"forestgreen\", lwd=2);\n";
	print TMPRSCRIPT "\tlines(x=position, y=atWeightedMedian, col=\"royalblue4\", lwd=2);\n";
	print TMPRSCRIPT "\tabline(h=0);\n";
	print TMPRSCRIPT "\tlegend(\"topright\", legend, cex=1, bty=\"n\", lty=c(1, 1), lwd=c(2, 2), col=legendColor);\n";
	print TMPRSCRIPT "dev.off()\n";
	close (TMPRSCRIPT);
	print "DONE\n";
	$cmd="R --vanilla --args < $prefix.temp.R";
	print "--> Running R script for plotting: $cmd ..";
	system ($cmd) == 0 || die "Unable to run R script: error code $?\n$!\n";
	print "DONE\n";
}


