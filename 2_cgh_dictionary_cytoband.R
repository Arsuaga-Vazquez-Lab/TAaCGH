#################################################
# USAGE
#
# Given a dataset, this code will give the 0-based start and stop indices for
# the sections in the chromosome arms.
#
# Input: aCGH file with no missing values formatted as specified in README.txt
# If the name of the dataset is "set", the name of the file should end with "data_full.txt"
# that is set_data_full.txt
#
# Output: Tab delimited file named  set_dict_cyto.txt with the following columns:
# "Chrom", "Arm", "Beg", "End", "Length", "Segment", "bpStart", "bpEnd", "CytoStart", "CytoEnd"
# where "bpStart", "bpEnd" are the start-index and stop-index repectively for the section.
# The index is 0-based, so downstream analysis with python or other 0-start languages
# can use the dictionary as is. Otherwise, in 1-start languages (R, etc.), the code will
# need to account for this.

#################################################
# COMMAND LINE ARGUMENTS
# 4. dataSet: short name for dataSet (e.g. set)
# 5. numParts: number of parts to split the dictionary. Usually 8
# 6. action: arms, sections
# 7. segLength: Section size (Best: 20 to 50). This will be the minimum number of probes
#       needed to run a specific arm. In the case of arms it will take the full arm but 
#       a value is needed (greater than the smaller arm number of clones)
# 8. subdir  a directory within /dataSet dir to read the dictionaries

# EXAMPLE
# Note: use vanilla when testing and slave when ready to avoid the log
# R --vanilla --args set 8 arms 20 arms < 2_cgh_dictionary_cytoband.R

library(gtools)

# Get the command line arguments
args = commandArgs();

dataSet <- args[4];
numParts <- as.integer(args[5]);
action <- args[6];
segLength <- as.integer(args[7]);
subdir <- args[8];

# to run locally
#dataSet <- "horlings";
#numParts <- 7;
#action <- "arms";
#segLength <- 20;
#subdir <- "arms";

###############################
# FUNCTIONS NEEDED

slice <- function(x,n) {
	N <- length(x);
	lapply(seq(1, N, n), function(i) x[i:min(i+n-1,N)]);
}


begPath <- "~/Research";

###############################
# READ aCGH FILE

dataFile <- paste(dataSet, "data", "full.txt", sep="_");
dataPath <- paste(begPath, "Data", dataSet, dataFile, sep="/");
# path to create a duplicate
dupFile  <- paste(dataSet, "data", "orig.txt", sep="_");
dupPath <- paste(begPath, "Data", dataSet, dupFile, sep="/");
print(dataPath);

data <- read.table(dataPath, header=T, sep='\t', comment.char="", stringsAsFactors = FALSE);

# Making sure the dataset was in the right order. Will replace original file
# saving a copy of the original file under set_data_orig.txt
# if Chrom has characters (like "X") will read as character variable and 
# will have produce a weird order (I believe it doesn't matter)
#TODO will be nice to have a proper routine to check for format complience
write.table(data,dupPath,sep='\t',row.names = FALSE)
data<-data[order(data$Chrom, data$Arm, data$bp),]
write.table(data,dataPath,sep='\t',row.names = FALSE)

###############################
# BEGIN PROGRAM


# Strip out the full list of chromosomes and arms
chromList <- data$Chrom;
armList <- data$Arm;

# Determine the unique chromosomes
chroms <- unique(chromList);
arms <- c("p", "q");

dict <- c();

for(chrom in chroms) {
	for(arm in arms) {
		# Collect the indices
		indices <- intersect(which(data$Chrom == chrom), which(data$Arm == arm));
		
		if(length(indices) != 0) {
			# The -1s make it so the indices are 0-based for later use
			# in the python scripts hom_stats.py and net_stats.py
			beg <- indices[1] - 1;
			end <- indices[length(indices)] - 1;
			armLength <- end - beg;
			#if(armLength > minProbes) {
            if(armLength > segLength) {

				if(action == "arms") {
					dict <- rbind(dict, c(chrom, arm, beg, end, end-beg, 1));
				} else if(action == "segments") {
					# Initialize variables to track location of segment for while loop
					segNum <- 1;
					#segLength <- 20;
					segAdvance <- round(segLength/2);
					segBegIndex <- beg;
					segEndIndex <- beg + segLength - 1;
					
					while(segEndIndex <= end) {
						# If there is a buffer from the end bigger than 15, add a non-last segment
						#if(abs(end - segEndIndex) > 9) {
                        if(abs(end - segEndIndex) > segAdvance-1) {
                            # These indices are getting things from R array so need to be base 1.
							bpStart <- data$bp[segBegIndex+1];
							bpEnd <- data$bp[segEndIndex+1];
							cytoStart <- toString(data$Cytoband[segBegIndex+1]);
							cytoEnd <- toString(data$Cytoband[segEndIndex+1]);
							
							dict <- rbind(dict, c(chrom, arm, segBegIndex, segEndIndex, segEndIndex-segBegIndex, segNum, bpStart, bpEnd, cytoStart, cytoEnd));
						#	dict <- rbind(dict, c(chrom, arm, segBegIndex, segEndIndex, segEndIndex-segBegIndex, segNum));
							
							# Update stuff
							segNum <- segNum + 1;
							segBegIndex <- segBegIndex + segAdvance;
							segEndIndex <- segBegIndex + segLength - 1;
						} else {
							# These indices are getting things from R array so need to be base 1.
							bpStart <- data$bp[segBegIndex+1];
							bpEnd <- data$bp[end+1];
							cytoStart <- toString(data$Cytoband[segBegIndex+1]);
							cytoEnd <- toString(data$Cytoband[end+1]);
							
							# Otherwise, this should be the last segment
							dict <- rbind(dict, c(chrom, arm, segBegIndex, end, end-segBegIndex, segNum, bpStart, bpEnd, cytoStart, cytoEnd));
						#	dict <- rbind(dict, c(chrom, arm, segBegIndex, end, end-segBegIndex, segNum));
							
							segBegIndex <- segBegIndex + segAdvance;
							segEndIndex <- segEndIndex + segLength - 1;
						}
					}
        }
			}
		}
	}
}

if(action == "arms") {
	# Change column names
	colnames(dict) <- c("Chrom", "Arm", "Beg", "End", "Length","OneSeg");
} else if(action == "segments") {
	# Change column names
	colnames(dict) <- c("Chrom", "Arm", "Beg", "End", "Length", "Segment", "bpStart", "bpEnd", "CytoStart", "CytoEnd");
#	colnames(dict) <- c("Chrom", "Arm", "Beg", "End", "Length", "Segment");
}

# Split up the dictionary
indices <- c(1:nrow(dict));
partsList <- slice(indices, ceiling(nrow(dict) / numParts));
for(i in c(1:length(partsList)))
{
	partDictFile <- paste(dataSet, '_dict_', i, '.txt', sep='');
	partDictPath <- paste(begPath, 'Data', dataSet, partDictFile, sep='/');
	
	print(dict[partsList[[i]], ]);
	
	write.table(dict[partsList[[i]], ], partDictPath, sep='\t', row.names=FALSE);
}

# Setup and write the dictionary file
dictFile <- paste(dataSet, "dict_cyto.txt", sep="_");
dictPath <- paste(begPath, "Data", dataSet, subdir, dictFile, sep="/");

write.table(dict, dictPath, sep='\t', row.names=FALSE);