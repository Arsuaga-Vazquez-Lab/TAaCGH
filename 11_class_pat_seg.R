# TODO: need to create a warning when saving phenotype file to avoid overwritting variables
# the script compares the homology curve from a patient against the average homology
# curve for test and for control using the leave-one-out procedure
# As a result, the phenotype file gets a new variable with every patient classified 
# for the section used as input

# OUTPUT
# will create additional variable in the phenotype file. The name of the variables
# will use input parameters (for example HER2_17qseg1_sig
# the values are: 1 if the segment was found significant for the patient, 0 otherwise, NA


param <- "B0";
phenotype <- "ER_pos";
dataSet <- "climent";
subdir <- "17qs1tos4"
action<- "sect";
arm="5p"
seg="seg3"
outliers <- "no";


###############################
# READ FILES

# Establish the beginning path
begPath <- "~/Research";

# The source code for the permutations
srcPath <- paste(begPath, "Code", "functions_sig.R", sep="/");
source(srcPath);

# Read the phenotype data
phenFile <- paste(dataSet, "phen.txt", sep="_");
phenPath <- paste(begPath, "Data", dataSet, phenFile, sep="/");
phenData <- read.table(phenPath, header=TRUE, sep="\t");
print(paste("phenPath",phenPath, sep=" "));

# Read data with homology
dataFile <- paste(param,"_2D_",dataSet,"_",action,"_",arm,"_", seg, ".txt", sep="");
#dataPath <- paste(begPath, "Results", dataSet, action, "2D/Homology",dataFile, sep="/");
dataPath <- paste(begPath, "Results", dataSet, subdir, "2D/Homology",dataFile, sep="/");
print(paste("dataPath homology",dataPath, sep=" "));

# Need to determine the maximum number of columns in the file ahead of time
# because the file is jagged. Not squared
nColData <- max(count.fields(dataPath, sep="\t"));

# Read in parameter data
data <- read.table(dataPath, header=FALSE, sep="\t", fill=TRUE, col.names=1:nColData);

# Fill in the NAs in the jagged array with a parameter appropriate value.
# For B0 and CC that value is 1 because that is the minimum/maximum, respectively.
# For B1 that value is 0.
# For D that value is the end index - start index in chrDict (0-based from cgh_dictionary.py).
 if (identical(param,"CC") || identical(param,"B0")) {
  fillValue <- 1.0;
  data[is.na(data)] <- fillValue;
} else if (identical(param,"B1")) {
  fillValue <- 0.0;
  data[is.na(data)] <- fillValue;
}

# getting the indices for each phenotype group 
phen1indices <- which(phenData[,phenotype] == 1);
phen2indices <- which(phenData[,phenotype] == 0);
phen1num<-length(phen1indices);
phen2num<-length(phen2indices);

print(paste("phenotype", phenotype, sep=" "));
print(paste("There are ", phen1num, " phenotype1", sep=""));
print(paste("There are ", phen2num, " phenotype2", sep=""));

# Initializing variables
phen1indices_noOut<-phen1indices;
phen2indices_noOut<-phen2indices;
phen1num_noOut<-phen1num;
phen2num_noOut<-phen2num;
outlier<-paste("out_",arm, sep="");

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

# Initializing the new variable which will indicate if the patient is significant for the segment
phenData$newVar <- NA
sumsExp <- colSums(exp)
AvgExp <- sumsExp / phen1num_noOut
sumsCon <- colSums(con)
AvgCon <- sumsCon / phen2num_noOut

# sum of squares leave-one-out
# first for the patients in the test group
for (i in phen1indices_noOut) {
  AvgExpLO <- (sumsExp-data[i,]) / (phen1num_noOut-1)
  sumSqVsExp <- sum((AvgExpLO - data[i,])^2)
#  phenData$sumSqVsExp[i] <- sum((AvgExpLO - data[i,])^2)
  sumSqVsCon <- sum((AvgCon - data[i,])^2)
#  phenData$sumSqVsCon[i]<- sum((AvgCon - data[i,])^2)
  if (sumSqVsExp > sumSqVsCon) {
#  if (phenData$sumSqVsExp[i] > phenData$sumSqVsCon[i]) {
    phenData$newVar[i] <- 0
  } else {
    phenData$newVar[i]<- 1
  }
} 
# now for the patients in the control group
for (i in phen2indices_noOut) {
  sumSqVsExp <- sum((AvgExp - data[i,])^2)
#  phenData$sumSqVsExp[i] <- sum((AvgExp - data[i,])^2)
  AvgConLO <- (sumsCon -data[i,]) / (phen2num_noOut-1)
  sumSqVsCon <- sum((AvgConLO - data[i,])^2)
#  phenData$sumSqVsCon[i]<- sum((AvgConLO - data[i,])^2)
  if (sumSqVsExp > sumSqVsCon) {
#  if (phenData$sumSqVsExp[i] > phenData$sumSqVsCon[i]) {
    phenData$newVar[i] <- 0
  } else {
    phenData$newVar[i]<- 1
  }
} 

test_stat <- sum(phenData$newVar[phen1indices_noOut])/phen1num_noOut
print(paste( phenotype, "_", arm ,seg, '_sig_perc =', test_stat, sep=""))

ctrl_stat <- sum(phenData$newVar[phen2indices_noOut])/phen2num_noOut
print(paste("no", phenotype, "_", arm ,seg, '_sig_perc =', ctrl_stat, sep=""))

cols <- ncol(phenData)
colnames(phenData)[cols] <- paste(phenotype, "_", arm, seg, "_sig", sep="")

# Update phen file
write.table(phenData, file=phenPath, row.names=F, sep='\t')
