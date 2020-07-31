warings    <- warnings();
args       <- commandArgs(TRUE);
name       <- args[1];
out        <- args[2];
mainTitle  <- args[3];
span       <- args[4]; # encode or genomewide
retBinInfo <- args[5]; # 0 = draw plots, 1 = return bins and counts 

if(span == "encode"){
	color <- "blue";
}else{
	color <- "black";
}
########## USER DEFINED FUNCTIONS #################
# Draw Density line graphs
draw_density_plots <- function(difference) {
	# # Create the png filename
	# histname <- noquote(strsplit(out,".png")[[1]]);
	# pngfile  <- paste(histname,"_density.png",sep='');
	# # Start PNG device driver to save output to figure.png
	# png(pngfile,height=3600,width=3600, res=300);

	#Draw plots
	plot (density(difference, from =0, to=0.5), xlab="(Midpoint of TAD - Midpoint of regulatory element)/TAD size", col=color, main=mainTitle, lwd=2, ylim=c(ymin,ymax), cex.main=1);
	abline(v=0,col="orange",lwd=2);
	abline(v=0.5 , col = "gray", lty = "dotted");
	legend("bottomright",span,col=color, fill=color, bty="l", ncol=1, cex=0.85);
}

# Draw histogram graphs
draw_histogram <- function(difference) {
	#Draw plots
	## its plotting till 0.49 as that is the mid of the last bin.
	allhist       <- hist(difference, right=FALSE, plot=FALSE);
	allhistmids   <- allhist$mids;
	allhistbreaks <- allhist$breaks;
	allhistcounts <- append(allhist$counts,0, after=length(allhist$counts));

	# # Create the png filename
	# histname <- noquote(strsplit(out,".png")[[1]]);
	# pngfile  <- paste(histname,"_frequency.png",sep='');
	# # Start PNG device driver to save output to figure.png
	# png(pngfile,height=3600,width=3600, res=300);

	plot(allhistbreaks, allhistcounts, main=mainTitle, lwd=2, cex.main=1, type="l",xlab="(Midpoint of TAD - Midpoint of regulatory element)/TAD size",ylim=c(0, max(allhistcounts)+max(allhistcounts)/10), xlim=c(0,0.49), col=color);
	abline(v=0,col="orange",lwd=2);
	abline(h=max(allhistcounts),col="gray", lwd=0.5, lty=3 );
	allymax  <- allhistbreaks[allhistcounts==max(allhistcounts)];
	text(allymax,max(allhistcounts), max(allhistcounts), col=2, adj=c(-.1,-.1),cex=0.5);
	legend("top",span, bty="l", ncol=1, col=color, fill=color, cex=0.75);
}

# Draw genomic distance graphs
draw_distance_plots <- function(difference, xlabel) {
	#Draw plots
	## its plotting till 0.49 as that is the mid of the last bin.
	allhist       <- hist(difference, right=FALSE, plot=FALSE);
	allhistmids   <- allhist$mids;
	allhistbreaks <- allhist$breaks;
	allhistcounts <- append(allhist$counts,0, after=length(allhist$counts));

	plot(allhistbreaks, allhistcounts, main=mainTitle, lwd=2, cex.main=1, type="l",xlab=xlabel, ylim=c(0, max(allhistcounts)+max(allhistcounts)/10), col=color);
	abline(v=0,col="orange",lwd=2);
	abline(h=max(allhistcounts),col="gray", lwd=0.5, lty=3 );
	allymax  <- allhistbreaks[allhistcounts==max(allhistcounts)];
	text(allymax,max(allhistcounts), max(allhistcounts), col=2, adj=c(-.1,-.1),cex=0.5);
	legend("top",span, bty="l", ncol=1, col=color, fill=color, cex=0.75);
}

###########################################################
# Read values from tab-delimited input file
# F1 = REVERSE Fragment
# F2 = FORWARD Fragment
# R1 = TADs
             
					   # __
						# /  \
					 # /    \
					# /      \
			   # /        \
		  	# /          \
			 # /     i      \
			# /     / \      \
		 # /     /   \      \
		# /     /     \      \
   # /     /       \      \
  # /     /         \      \
 # /     /           \      \
# +-----|-------|-----|------+
# Ts    Rm      Tm    Fm     Te


# getting the file name
inputfile <-paste(name,sep='');
# Read values from tab-delimited input file
# chr1    770137  1250137 chr1    794032  794124  WE      1000    .       794032  794124  255,252,4       92
# chr1    770137  1250137 chr1    800800  801452  WE      1000    .       800800  801452  255,252,4       652
# chr1    770137  1250137 chr1    801544  801849  WE      1000    .       801544  801849  255,252,4       305
# chr1    770137  1250137 chr1    823384  823549  WE      1000    .       823384  823549  255,252,4       165
# chrX    154106806       154946806       chrX    154697946       154697947       ENST00000447347.1       0       +       1
# chrX    154106806       154946806       chrX    154718983       154718984       ENST00000452506.1       0       +       1
# chrX    154106806       154946806       chrX    154722149       154722150       ENST00000449645.2       0       -       1
# chrX    154106806       154946806       chrX    154774937       154774938       ENST00000487422.1       0       -       1

# Get the genome-wide or encode data
alldata <-read.table(inputfile,header=FALSE,sep="\t", na.strings="NA");
alltadChrom <- alldata$V1;
alltadStart <- alldata$V2;
alltadEnd   <- alldata$V3;
alltssChrom <- alldata$V4;
alltssStart <- alldata$V5;
alltssEnd   <- alldata$V6;
alltssName  <- alldata$V7;

# get tads stats
minTadSize <- min(alltadEnd - alltadStart);
maxTadSize <- max(alltadEnd - alltadStart);
avgTadSize <- round(mean(alltadEnd - alltadStart));
cat("\nminTadSize=", minTadSize,"\nmaxTadSize=", maxTadSize, "\navgTadSize=", avgTadSize, "\n\n");

# handle strand column
# if its completely empty 
if(is.null(alldata$V9)){ 
	alltssStrand <- "+";
}else{
	alltssStrand<- alldata$V9;
}

# if there is no strand
if(alltssStrand == '.' || alltssStrand == ''){
	alltssStrand <- '+';
}

# Calculate the difference for the looping data
allposTssStart <- alltssStart[alltssStrand=='+'];
allposTssEnd   <- alltssEnd[alltssStrand=='+'];
allnegTssStart <- alltssEnd[alltssStrand=='-'];
allnegTssEnd   <- alltssStart[alltssStrand=='-'];

allposTadStart <- alltadStart[alltssStrand=='+'];
allposTadEnd   <- alltadEnd[alltssStrand=='+'];
allnegTadStart <- alltadStart[alltssStrand=='-'];
allnegTadEnd   <- alltadEnd[alltssStrand=='-'];

allposTadMid <- allposTadStart + round((allposTadEnd - allposTadStart)/2);
allnegTadMid <- allnegTadStart + round((allnegTadEnd - allnegTadStart)/2);

allposTssMid <- allposTssStart + round((allposTssEnd - allposTssStart)/2);
allnegTssMid <- allnegTssStart + round((allnegTssEnd - allnegTssStart)/2);



## Calculate the distance from the mid point
allposMidDiff <- abs(allposTadMid - allposTssMid);
allnegMidDiff <- abs(allnegTadMid - allnegTssMid);
allMiddifference <- c(allposMidDiff, allnegMidDiff);

## Calculate the distance from the start or end which ever is small
allposStartDiff <- abs(allposTadStart - allposTssMid);
allnegStartDiff <- abs(allnegTadStart - allnegTssMid);
allposEndDiff   <- abs(allposTadEnd   - allposTssMid);
allnegEndDiff   <- abs(allnegTadEnd   - allnegTssMid);


allPosPosBoundaryDiff <- allposStartDiff[allposStartDiff <= allposEndDiff];
allPosNegBoundaryDiff <- allposEndDiff[allposStartDiff > allposEndDiff];

allNegPosBoundaryDiff <- allnegStartDiff[allnegStartDiff <= allnegEndDiff];
allNegNegBoundaryDiff <- allnegEndDiff[allnegStartDiff > allnegEndDiff];

allBoundarydifference <- c(allPosPosBoundaryDiff, allPosNegBoundaryDiff, allNegPosBoundaryDiff, allNegNegBoundaryDiff);

# Single Sided graphs
allposDiff <- abs(allposTadMid - allposTssMid)/(allposTadEnd - allposTadStart);
allnegDiff <- abs(allnegTadMid - allnegTssMid)/(allnegTadEnd - allnegTadStart);

# Get the difference
alldifference <- c(allposDiff, allnegDiff);

# Make all the values greater than 0.5 to 0.5 and less than 0 to 0
alldifference[alldifference > 0.5] <- 0.5;
alldifference[alldifference < 0]   <- 0;

allyd    <- density(alldifference, from =0, to=0.5);
allydmin <- min(allyd$y);
allydmax <- max(allyd$y);

ymin <- min(allydmin);
ymax <- max(allydmax);

# One sided Graphs
if(retBinInfo == 0){
	# Create the png filename
	histname <- noquote(strsplit(out,".png")[[1]]);
	pngfile  <- paste(histname,"_distance_freq_density.png",sep='');
	# Start PNG device driver to save output to figure.png
	png(pngfile,height=3600,width=3600, res=300);
	par(mfrow=c(2,2)); # 2 row 2 plots
	draw_distance_plots(allBoundarydifference, "(Distance of midpoint of regulatory element from the closest TAD boundary)");
	draw_distance_plots(allMiddifference, "(Distance of midpoint of regulatory element from the midpoint of the TAD)");
	draw_density_plots(alldifference);
	draw_histogram(alldifference);
} else if(retBinInfo == 1){
	bins <- seq(0,0.525, by=0.025);
	allhist       <- hist(alldifference, right=FALSE, plot=FALSE, breaks=bins);
	allhistbreaks <- allhist$breaks[1:(length(allhist$breaks) -2)]; # remove the last bin
	allhistcounts <- allhist$counts[-length(allhist$counts)]; # counts are always one less than the bins. 

	#Convert counts into percentages.
	allhistpercs <- allhistcounts*100/sum(allhistcounts);
	cat(allhistbreaks,"|",allhistpercs);
}


