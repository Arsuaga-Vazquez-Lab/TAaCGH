# This program is "new" because it allows you to work with full arms as well as sections
# and is out becuase allows you to handle outliers
#################################################
# USAGE
# 
# Given the jagged array files created by 4_hom_stats_parts.py,
# this code will compute a p-value for each chromosome arm and dimension under investigation.
# We use permutation tests to arrive at the p-value, and the permutation aspects are located
# in functions_sig.R.
#
# Note: The result is the collection of uncorrected p-values, FDR will be used in scripts to follow
# Make sure to have adifferent name for the dataset in case you want to runit by arms and sections
# for example: SET_arms and SET_sect as otherwise it might overwrite the first output
#
# Input: The following tab delimited text files.
# 1. The dictionary files from 2_cgh_dictionary_cytoband.R.
# 2. The phenotype file where the patient names matches the order from the aCGH file with the log2 ratios.
# A column with the phenotype (ERBB2, basal, test, sim, etc) coded with 0's and 1's. The column name will be the
# argument "phenotype"
# 3. The jagged array (from 4_hom_stats_parts.py) containing the value of the parameter (B0, B1) for all the patient
# simplicial complexes/networks over the epsilon filtrations. The files should be under 
# ~/Results/SET/2D/Homology/B0_2D_SET_chrArm_seg.txt for a specific aCGH file named SET and a specific chromosome-arm
# combination "chrArm"
#
# Output:
# A tab delimited text file with 12 columns (# "Chr", "Arm", "Dimension", "Segment",
# "Idx0.Beg", "Idx0.End", "bpStart", "bpEnd", "CytoStart", "CytoEnd","P.Value", "Test.Stat") which
# gives the p-value and test statistic for each chromosome arm for each computed dimension. 
# "Idx0.Beg" and "Idx0.End" are the beginning and ending positions in the CGH file for the sections (0 start)
# Output will be saved under ~/Research/Results/SET/significance/pvals

#################################################
# COMMAND LINE ARGUMENTS
# 4. param (B0, B1, D, CC)
# 5. phenotype (ERBB2, basal, test, sim, etc)
# 6. dataSet (SET, etc)
# 7. partNum (an integer with the part from the dictionary)
# 8. action: arms, sections
# 9. outliers to remove: yes, no

# EXAMPLE
# R --vanilla --args B1 Basal horlings_arms 3 arms yes< 5_sig_pcalc_parts_new_out.R
# R --vanilla --args B0 Luminal_A bergamaschi1pMADMA3 1 sections no< 5_sig_pcalc_parts.R
#  

# Get the command line arguments
args = commandArgs();

param <- args[4];
phenotype <- args[5];
dataSet <- args[6];
partNum <- args[7];
action <- args[8];
outliers <- args[9];
 
# for debugging purposes only
# param <- "B0";
# phenotype <- "TP53_mut";
# dataSet <- "bergamaschi_sect";
# partNum <- "11";
# action<- "arms";
# outliers <- "yes";

###############################
# READ FILES

# Establish the beginning path
begPath <- "~/Research";

# Determine "Homology" or "Networks" based on param
if (identical(param,"D") || identical(param,"CC")) {
	paramType <- "Networks";
	shortParamType <- "net";
} else {
	paramType <- "Homology";
	shortParamType <- "hom"
}

# The source code for the permutations
srcPath <- paste(begPath, "Code", "functions_sig.R", sep="/");
source(srcPath);

# Read the chromosome dictionary data
chrDictFile <- paste(dataSet, "_dict_", partNum, ".txt", sep="");
chrDictPath <- paste(begPath, "Data", dataSet, chrDictFile, sep="/");
chrDict <- read.table(chrDictPath, header=TRUE, sep="\t");
print(chrDict);

# Read the phenotype data
phenFile <- paste(dataSet, "phen.txt", sep="_");
phenPath <- paste(begPath, "Data", dataSet, phenFile, sep="/");
phenData <- read.table(phenPath, header=TRUE, sep="\t");

###############################
# BEGIN PROGRAM

# getting the indices for each phenotype group 
phen1indices <- which(phenData[,phenotype] == 1);
phen2indices <- which(phenData[,phenotype] == 0);
phen1num<-length(phen1indices);
phen2num<-length(phen2indices);

print(paste("phenotype", phenotype, sep=" "));
print("phenData$phenotype", phenData[,phenotype]);

print(paste("There are ", phen1num, " phenotype1", sep=""));
print(paste("There are ", phen2num, " phenotype2", sep=""));


# Write path for p-values and header
uncorrFile <- paste(param, "_", phenotype, "_", dataSet, "_pvals_", partNum, ".txt", sep="");
uncorrFolder <- paste(begPath, "Results", dataSet, "significance", "pvals", sep="/");

if(!file.exists(uncorrFolder)) {
	dir.create(uncorrFolder, recursive=TRUE);
}

uncorrPath <- paste(uncorrFolder, uncorrFile, sep="/");

# Create a container for the p-value list
pvals <- c();

for(i in c(1:nrow(chrDict)))
{
	chr <- chrDict[i,1];
	arm <- as.vector(chrDict[i,2]);
	chrArm <- paste(chr, arm, sep="");
	seg <- chrDict[i,6];
	beg <- chrDict[i,3];
	end <- chrDict[i,4];
    if (action=="sections") {
		bpStart <- chrDict[i,7];
		bpEnd <- chrDict[i,8];
		CytoStart <- as.vector(chrDict[i,9]);
		CytoEnd <- as.vector(chrDict[i,10]);
    }
	print(paste("On chromosome ", chrArm, " segment ", seg, sep=""));
	
	dim<-"2";

	print(paste("On dimension ", dim, sep=""));
		
	# Read in data regardless of phenotype
	dataFile <- paste(param, "_", dim, "D_", dataSet, "_", chrArm, "_seg", seg,".txt", sep="");
	dataPath <- paste(begPath, "/Results/", dataSet, "/", dim, "D/", paramType, "/", dataFile, sep="");
		
	print(dataPath);
		
	# Need to determine the maximum number of columns in the file ahead of time
	# because, again, R is strange about reading in jagged arrays.
	nColData <- max(count.fields(dataPath, sep="\t"));
		
	# Read in parameter data
	data <- read.table(dataPath, header=FALSE, sep="\t", fill=TRUE, col.names=1:nColData);
		
	# Fill in the NAs in the jagged array with a parameter appropriate value.
	# For B0 and CC that value is 1 because that is the minimum/maximum, respectively.
	# For B1 that value is 0.
	# For D that value is the end index - start index in chrDict (0-based from cgh_dictionary.py).
	if (identical(param,"D")) {
		fillValue <- chrDict[i,4] - chrDict[i,3]
		data[is.na(data)] <- fillValue;
	} else if (identical(param,"CC") || identical(param,"B0")) {
		fillValue <- 1.0;
		data[is.na(data)] <- fillValue;
	} else if (identical(param,"B1")) {
		fillValue <- 0.0;
		data[is.na(data)] <- fillValue;
	}
    
    
    # Initializing variables
    phen1indices_noOut<-phen1indices;
    phen2indices_noOut<-phen2indices;
    phen1num_noOut<-phen1num;
    phen2num_noOut<-phen2num;
    outlier<-paste("out_",chrArm, sep="");
    
    # Remove outliers from arm
	if (outliers=="yes") {
	    if ( (outlier %in% colnames(phenData))==TRUE) {
    	    outlierIndx<-which(phenData[,outlier]==1);
        	phen1indices_noOut<-setdiff(phen1indices, outlierIndx);
	        phen1num_noOut<-length(phen1indices_noOut);
    	    phen2indices_noOut<-setdiff(phen2indices, outlierIndx);
        	phen2num_noOut<-length(phen2indices_noOut);
	        print(paste("After removing outliers in ",outlier,":", sep=""));
    	    print(paste("There are ", phen1num_noOut, " patients for test", sep=""));
        	print(paste("There are ", phen2num_noOut, " patients for control", sep=""));
    	}
    }

	# Separate by phenotype
	exp <- data[phen1indices_noOut,];
	con <- data[phen2indices_noOut,];

	# Compute the pvalue
	rndData <- randomize(exp, con, fillValue);
	pval <- rndData[[1]];
	stat <- rndData[[2]];
	sumPos <- rndData[[3]];
	sumNeg <- rndData[[4]];
		
	print(pval);
	print(stat);
	print(partNum);
		
	# Put the pvals in the container
    if (action=="sections") {
		pvals <- rbind(pvals, c(chr, arm, dim, seg, phen1num_noOut, phen2num_noOut, beg, end, bpStart, bpEnd, CytoStart, CytoEnd, pval, stat, sumPos, sumNeg));
    } else {
        pvals <- rbind(pvals, c(chr, arm, dim, seg, phen1num_noOut, phen2num_noOut, beg, end, pval, stat, sumPos, sumNeg));
    }
	
}

# Rename columns and write the file
if (action=="sections") {
    colnames(pvals) <- c("Chr", "Arm", "Dimension", "Segment", "TestNum", "CtrlNum", "Idx0.Beg", "Idx0.End", "bpStart", "bpEnd", "CytoStart", "CytoEnd","P.Value", "Test.Stat", "AreaPos", "AreaNeg");
} else {
    colnames(pvals) <- c("Chr", "Arm", "Dimension", "Segment", "TestNum", "CtrlNum", "Idx0.Beg", "Idx0.End","P.Value", "Test.Stat", "AreaPos", "AreaNeg");
}
write.table(pvals, file=uncorrPath, row.names=F, sep='\t')
