# This programs creates the transposed version of the data set
# to feed into hom_stats.py.

# Input: xxxx_data_full.txt
# Output: xxxx_data.txt


#################################################
# COMMAND LINE ARGUMENTS
# 4. loc (local, gonzalez)
# 5. dataSet (sj, sim6, simC3 climent, etc.)


# EXAMPLE
# Note: use vanilla when testing and slave when ready to avoid the log
# R --slave --args simC2  < Transposed_aCGH.R
# R --vanilla --args set < 3_Transposed_aCGH.R

# Get the command line arguments
args = commandArgs();

dataSet <- args[4];



 CGH_start <- 6;
###############################
# READ FILES
begPath <- "~/Research";
	
	dataFile <- paste(dataSet, "data", "full.txt", sep="_");
	dataPath <- paste(begPath, "Data", dataSet, dataFile, sep="/");
	
	print(dataPath);
	
	data <- read.table(dataPath, header=T, sep='\t', comment.char="");
	
	###############################
	# BEGIN PROGRAM
	
	# Keeping Chrom and Arm info in a transposed data frame
	dataT <- as.data.frame(t(data$Chrom));
	dataT <- rbind(dataT, t(as.character(data$Arm)));
	
	# Getting the number of rows and columns from data_full
	cols <- ncol(data);
	rows <- nrow(data);
	
	# Adding to the data frame the transposed patients aCGH info
	dataT <- rbind(dataT, t(data[,c(CGH_start:cols)]));
	
	# Setup and write the dictionary file
	newFile <- paste(dataSet, "data.txt", sep="_");
	newPath <- paste(begPath, "Data", dataSet, newFile, sep="/");
	
	write.table(dataT, newPath, sep='\t', row.names=FALSE, col.names=FALSE);	
