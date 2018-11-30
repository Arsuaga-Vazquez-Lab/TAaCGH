# This script will classify a patient as significant for the center of mass 
# if the average for the arm (CM of the patient) lays outside de
# confidence interval from the control group of the phenotype

# The confidence interval is from a t distribution: 
# avg +/- t*sd/sqrt(n)   t=alfa/2 percentile from a t with n-1 d.f.

# INPUT
# size, mean and standard deviation for the control group from
# the file created by 9_mean_diff_perm_NoOut
# TODO: READ SIZE, MEAN AND SD FROM FILE INSTEAD AS ARGUMENT

# OUTPUT
# a new variable will be create in the phenotype file SET_phen.txt


# Here is the info for Horlings ER+, n_ctrl=28
#  16p mean_ctrl=0.02, sd=0.08
#  16q mean_ctrl= 0, sd=0.11
  
# Only for debugging purposes
 dataSet <- "horlings";
 phenotype <- "ER_pos";
 sig <- 0.05;
 chrom <- 16
 arm <- "q"
 n_ctrl <- 28
 mean_ctrl <- 0
 sd_ctrl <- 0.11
 # type: gain or del
 type <- "del"

# Working directory
begPath <- "~/Research";

# The following CGH_start is minus 1
CGH_start <- 6;


### INPUT/OUTPUT FILES####

# Read the phenotype data
phenFile <- paste(dataSet, "phen.txt", sep="_");
phenPath <- paste(begPath, "Data", dataSet, phenFile, sep="/");
phenData <- read.table(phenPath, header=TRUE, sep="\t");

# Get the CGH data
dataPath <- paste(begPath, '/Data/', dataSet, '/', dataSet, '_data_full.txt', sep='');
data <- read.table(dataPath, header=T, sep='\t', comment.char='');

######### BEGIN

# Subset and get only the data for the arm
dataIndices <- intersect(which(data$Chrom == chrom), which(data$Arm == arm));
arm_data <- data[dataIndices,CGH_start:ncol(data)];

# CM: center of mass by patient
CM_i <- colMeans(arm_data)

phenData$class <- 0

if (type=="gain") {
  upper <- mean_ctrl + qt(sig,df=n_ctrl-1, lower.tail=FALSE)*sd_ctrl/sqrt(n_ctrl)
  phenData$class[CM_i>upper] <- 1
}else if (type=="del") {
  lower <- mean_ctrl - qt(sig,df=n_ctrl-1)*sd_ctrl/sqrt(n_ctrl)
  phenData$class[CM_i<lower] <- 1
} else { print("incorrect input for argument type {gain, del}")}

colnames(phenData)[ncol(phenData)] <- paste(phenotype,"_CM_",chrom,arm,"_sig", sep="")

# Update phen file
write.table(phenData, file=phenPath, row.names=F, sep='\t')

