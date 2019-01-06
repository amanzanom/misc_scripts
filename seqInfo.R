

#===============================================================================
#   Author: Alejandro Manzano Marin
#
#   File: seqInfo.R
#   Date: 15-11-2012
#   Version: 1.0
#
#   Usage:
#      R --vanilla --args seqInfoOut.aux seqInfoPrefix [flx|plus] < pathTo/seqInfo.R
#
#    Description: This program will take the ouput *.aux file of seqInfo.pl and do plotting of the
#                 parameters chosen to be analyzed.
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




#Capture all data to the variables in R through reading the external file as a dataframe.
FILE<-as.character(commandArgs(trailingOnly=T)[1]);
outFile<-as.character(commandArgs(trailingOnly=T)[2]);
reads<-as.character(commandArgs(trailingOnly=T)[3]);
FILE;
nucleotideSeq<- read.table(FILE, header=FALSE, sep="\t", col.names=c("id", "length", "gc", "gcFirst", "gcSecond", "gcFplusS", "gcThird"), colClasses=c("character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"), na.strings="NULL");

pdf (paste(outFile, ".pdf", sep=""))
	
	#### Length hist
	x<-nucleotideSeq$length;
	xMean<- round(mean(x), 2);
	xStdDev<- round(sd(x), 2);
	xMin<- min(x);
	xMax<- max(x);
	if (reads == "flx"){
		x300<-length(x[x<400]);
		x400<-length(x[x>=400 & x<=500]);
		x500<-length(x[x>500]);
		legend <- c(paste("min= ", xMin), paste("max= ", xMax), paste("mean= ", xMean), paste("stdDev= ", xStdDev), paste("#Seqs= ", length(x)), paste("<=400bp= ", x300), paste("400-500bp= ", x400), paste(">500bp= ", x500)); 
	}
	if (reads == "plus"){
		x300<-length(x[x<600]);
		x400<-length(x[x>=600 & x<=700]);
		x500<-length(x[x>700]);
		legend <- c(paste("min= ", xMin), paste("max= ", xMax), paste("mean= ", xMean), paste("stdDev= ", xStdDev), paste("#Seqs= ", length(x)), paste("<700bp= ", x300), paste("700-800bp= ", x400), paste(">800bp= ", x500)); 
	}
	if (reads == "0"){
		legend <- c(paste("min= ", xMin), paste("max= ", xMax), paste("mean= ", xMean), paste("stdDev= ", xStdDev), paste("#Seqs= ", length(x))); 
	}
	hist(x, breaks=seq( from=xMin, to=xMax, by=(xMax-xMin)/100 ), xlim=c(xMin, xMax), main="Length distribution", xlab="Length (bp)");
	abline(v=c(xMean), col="red", lty=2);
	legend("topright", legend, cex=1, bty="n");
	
	#### GC frequencies
	if (length(nucleotideSeq$gc[!is.na(nucleotideSeq$gc)])){
		x<-nucleotideSeq$gc;
		xMean<- round(mean(x), 4);
		xStdDev<- round(sd(x), 4);
		xMin<- min(x);
		xMax<- max(x);
		legend <- c(paste("min= ", round(xMin, 4)), paste("max= ", round(xMax, 4)), paste("mean= ", xMean), paste("stdDev= ", xStdDev));
		hist(x, breaks=seq( from=xMin, to=xMax, by=(xMax-xMin)/100 ), xlim=c(0, 1), main="G+C distribution", xlab="G+C content");
		abline(v=c(xMean), col="red", lty=2);
		legend("topright", legend, cex=1, bty="n");
	}
	
	#### GC first hist frequencies
	if (length(nucleotideSeq$gcFirst[!is.na(nucleotideSeq$gcFirst)])){
		x<-nucleotideSeq$gcFirst;
		
		xMean<- round(mean(x), 4);
		xStdDev<- round(sd(x), 4);
		xMin<- min(x);
		xMax<- max(x);
		legend <- c(paste("min= ", round(xMin, 4)), paste("max= ", round(xMax, 4)), paste("mean= ", xMean), paste("stdDev= ", xStdDev));
		hist(x, breaks=seq( from=xMin, to=xMax, by=(xMax-xMin)/100 ), xlim=c(xMin, xMax), main="G+C 1st pos distribution", xlab="GC");
		abline(v=c(xMean), col="red", lty=2);
		legend("topright", legend, cex=1, bty="n");
	}
	
	#### GC second hist frequecies
	if (length(nucleotideSeq$gcSecond[!is.na(nucleotideSeq$gcSecond)])){
		x<-nucleotideSeq$gcSecond;
		xMean<- round(mean(x), 4);
		xStdDev<- round(sd(x), 4);
		xMin<- min(x);
		xMax<- max(x);
		legend <- c(paste("min= ", round(xMin, 4)), paste("max= ", round(xMax, 4)), paste("mean= ", xMean), paste("stdDev= ", xStdDev));
		hist(x, breaks=seq( from=xMin, to=xMax, by=(xMax-xMin)/100 ), xlim=c(xMin, xMax), main="G+C 2nd pos distribution", xlab="GC");
		abline(v=c(xMean), col="red", lty=2);
		legend("topright", legend, cex=1, bty="n");
	}
	
	#### GC first+second hist frequencies
	if (length(nucleotideSeq$gcFplusS[!is.na(nucleotideSeq$gcFplusS)])){
		x<-nucleotideSeq$gcFplusS;
		xMean<- round(mean(x), 4);
		xStdDev<- round(sd(x), 4);
		xMin<- min(x);
		xMax<- max(x);
		legend <- c(paste("min= ", round(xMin, 4)), paste("max= ", round(xMax, 4)), paste("mean= ", xMean), paste("stdDev= ", xStdDev)); 
		hist(x, breaks=seq( from=xMin, to=xMax, by=(xMax-xMin)/100 ), xlim=c(xMin, xMax), main="G+C 1st&2nd pos distribution", xlab="GC");
		abline(v=c(xMean), col="red", lty=2);
		legend("topright", legend, cex=1, bty="n");
	}
	
	#### GC third third hist frequencies
	if (length(nucleotideSeq$gcThird[!is.na(nucleotideSeq$gcThird)])){
		x<-nucleotideSeq$gcThird;
		xMean<- round(mean(x), 4);
		xStdDev<- round(sd(x), 4);
		xMin<- min(x);
		xMax<- max(x);
		legend <- c(paste("min= ", round(xMin, 4)), paste("max= ", round(xMax, 4)), paste("mean= ", xMean), paste("stdDev= ", xStdDev));
		hist(x, breaks=seq( from=xMin, to=xMax, by=(xMax-xMin)/100 ), xlim=c(xMin, xMax), main="G+C 3rd pos distribution", xlab="GC");
		abline(v=c(xMean), col="red", lty=2);
		legend("topright", legend, cex=1, bty="n");
	}
	
	#### G+C content vs Length
	if (length(nucleotideSeq$gc[!is.na(nucleotideSeq$gc)])){
		x<-nucleotideSeq$length;
		x<-cbind(x, nucleotideSeq$gc);
		plot(x[,1]/1000, x[,2], main="G+C content vs Length", xlab="Length (Kb)", ylab="G+C content", col=rgb(0,100,0,50,maxColorValue=255), pch=16, xlim=c(0,2200), ylim=c(0,1));
	}

dev.off();
