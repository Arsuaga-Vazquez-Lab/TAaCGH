# This program generates profiles by patient (ie: bp vs log2 ratio CGH) for an specific arm
# In one pdf you will have as many graphs as patients in an specific group
# Type in the data at the beginning with the chorm and arm for the graphs and phenotype
# Make sure you include all the patients you need in the group as this program will not read the phenotype
# modify the for loop for the name you used for the list of patients
# This program may also add 2 or 4 lvertical lines to show where the segments start, you have to input the position for them at the end 
# example in local: source("/Users/gina/Desktop/Gina/BioTec/Homology/Code Orig and Modified/myPaper/ind_prof_origpat_local.R");
begPath <- "~/Research/";
dataSet <- "horlings";
cghStart <- 6;
chrArm <- "4q";
chr <- "4";
arm <- "q";
phenotype <- "Basal";
subdir <- ""

dataPath <- paste(begPath, '/Data/', dataSet, '/', dataSet, '_data_full.txt', sep='');
data <- read.table(dataPath, header=T, sep='\t');
cgh <- data[, cghStart:ncol(data)];

numPats <- ncol(cgh);
# getting the indices for each phenotype group 

# Read the phenotype data
phenFile <- paste(dataSet, "phen.txt", sep="_");
phenPath <- paste(begPath, "Data", dataSet, phenFile, sep="/");
phenData <- read.table(phenPath, header=TRUE, sep="\t");
print(paste("phenPath",phenPath, sep=" "));
phen1indices <- which(phenData[,phenotype] == 1);
phen1num<-length(phen1indices);
print(paste("There are ", phen1num, " phenotype1", sep=""));



	caIndices <- intersect(which(data$Chrom == chr), which(data$Arm == arm));
	maxY <- max(cgh[caIndices,phen1indices])
	minY <- min(cgh[caIndices,phen1indices])
	
	# Write folder
	writeDir = paste(begPath, '/Data/', dataSet,'/',subdir, '/vis/ind_profiles_clouds/',phenotype,'/', sep="", collaspe=NULL);
	if(!file.exists(writeDir)) {
		dir.create(writeDir, recursive=TRUE);
	}
	
	print(paste("On ", chrArm, sep=''));

	profilePath <- paste(writeDir, phenotype,'_',chrArm,'.pdf', sep='');
	pdf(profilePath, width=9, height=6);
	par(mfrow=c(3,2), cex.main=1.5,cex.lab=1.2,mar=c(4,3,2,2) + 0.1);
# for only 2 plots
#	pdf(profilePath, width=9, height=3.3);
#	par(mfrow=c(1,2), cex.lab=1, cex.main=1.23, mar=c(4,3,2,2) + 0.1);

	for(p in phen1indices) {
		profile.title <- paste('Patient ', colnames(cgh[p]), ' on ', chrArm, sep='');
		profile.xlab <- 'bp';
		profile.ylab <- 'log2 ratio';
		
#		bp <- data[caIndices, 5]; #kwek has bp on the 5th column
		bp <- data[caIndices, 4]; # bp is the 4th column in horlings
		patData <- cgh[caIndices, p];
		
		plot(bp, patData, main=profile.title, xlab="", ylab="", ylim = c(minY, maxY), pch=20, type="o");
#       plot(bp, patData, main=profile.title, xlab="", ylab="", ylim = c(-10, 10), pch=20, type="o");
		abline(h=c(0), lty=2, col="black");
		# abline(v=72021300, lty=3, col="black");
		# abline(v=90704922, lty=3, col="black");
		# abline(v=81980978, lty=3, col="red");
		# abline(v=101292561, lty=3, col="red");
		
	}
	dev.off();


