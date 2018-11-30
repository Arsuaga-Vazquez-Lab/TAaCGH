# This program generates profiles by patient (ie: bp vs log2 ratio CGH) for a specific section provided inital and ending bp positions 
# In one pdf you will have as many graphs as patients in an specific group
# Type in the data at the beginning with the chorm and arm for the graphs and phenotype
# Make sure you include all the patients you need in the group as this program will not read the phenotype
# modify the for loop for the name you used for the list of patients
# This program may also add 2 or 4 vertical lines to show where the segments start, you have to input the position for them at the end
# example in local: source("/Users/gina/Desktop/Gina/BioTec/Homology/Code Orig and Modified/myPaper/ind_prof_origpat_local.R");

# TODO: Pass arguments in an efficient way

begPath <- "~/Research/";
dataSet <- "horlings";
cghStart <- 6;
chrArmSect <- "17q.s2";
chr <- "4";
arm <- "q";
starting <- 32489785  #starting position for the section in bp
ending <- 43339849  #ending position for the section in bp  
phenotype <- "HER2";
subdir <- "sect"

#for 17qs1
#starting <- 25440972
#ending <- 37812853
#for 17qs2
#starting <- 32489785
#ending <- 43339849

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



#caIndices <- intersect(which(data$Chrom == chr), which(data$Arm == arm));

	dataSect <- data[(data$Chrom==chr & data$bp>=starting & data$bp<=ending),]
  Mbp <- dataSect$bp/1000000
  cgh <- dataSect[,cghStart:ncol(data)]
	maxY <- max(cgh[,phen1indices])
	minY <- min(cgh[,phen1indices])
	
		
	# Write folder
	writeDir = paste(begPath, '/Data/', dataSet,'/',subdir, '/vis/ind_profiles/',phenotype,'/', sep="", collaspe=NULL);
	if(!file.exists(writeDir)) {
		dir.create(writeDir, recursive=TRUE);
	}
	
	print(paste("On ", chrArmSect, sep=''));

	profilePath <- paste(writeDir, phenotype,'_',chrArmSect,'.pdf', sep='');
	pdf(profilePath, width=9, height=6);
	par(mfrow=c(3,2), cex.main=1.5,cex.lab=1.2,mar=c(4,3,2,2) + 0.1);
# for only 2 plots
#	pdf(profilePath, width=9, height=3.3);
#	par(mfrow=c(1,2), cex.lab=1, cex.main=1.23, mar=c(4,3,2,2) + 0.1);

	for(p in phen1indices) {
		profile.title <- paste('Patient ', colnames(cgh[p]), ' on ', chrArmSect, sep='');
		profile.xlab <- 'Mbp';
		profile.ylab <- 'log2 ratio';
		
#		bp <- data[caIndices, 5]; #kwek has bp on the 5th column
#		bp <- data[caIndices, 4]; # bp is the 4th column in horlings
		patData <- cgh[, p];
		
		plot(Mbp, patData, main=profile.title, xlab="", ylab="", ylim = c(minY, maxY), pch=20, type="o");
#       plot(bp, patData, main=profile.title, xlab="", ylab="", ylim = c(-10, 10), pch=20, type="o");
		abline(h=c(0), lty=2, col="black");
		# abline(v=72021300, lty=3, col="black");
		# abline(v=90704922, lty=3, col="black");
		# abline(v=81980978, lty=3, col="red");
		# abline(v=101292561, lty=3, col="red");
		
	}
	dev.off();


