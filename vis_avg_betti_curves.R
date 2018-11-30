# Average B0 curves for Test and control compare with patient curve

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
phenotype <- "HER2";
dataSet <- "horlings";
subdir <- "17qs1tos4"
action<- "sect";
arm="17q"
seg="seg1"
outliers <- "no";
epsilon <- 0.05 # the increments used
patExp <- 18
patCon <- 37

# Epsilons used by case
# HER@ horlings 17qs1tos4 e=0.05

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

# chr <- sig_sections$Chr[i];
# arm <- as.vector(sig_sections$Arm[i]);
# chrArm <- paste(chr, arm, sep="");
# seg <- sig_sections$Segment[i];

print(paste("On chromosome ", arm, " segment ", seg,  sep=""));

# Plots will be saved under the following name
# Ensure the folder we are going to write to exists
curveFolder = paste(begPath, 'Results', dataSet, subdir, 'vis', 'class_curves', param, phenotype, sep='/');
if(!file.exists(curveFolder)) {
  dir.create(curveFolder, recursive=TRUE);
}
curveFile <- paste(param, '_', phenotype, '_2D_', dataSet, "_",action,'_', arm, '_', seg, '_patT_',patExp,'patC_',patCon,'.pdf', sep='');
curvePath <- paste(curveFolder, curveFile, sep='/');

# generating the ith curveMeans

#phen1curve <- curvesMeans_i$test #AvgExp
#phen2curve <- curvesMeans_i$control #AvgCon
yMax <- max(max(AvgExp), max(AvgCon));
xMax<-length(AvgExp)*epsilon-epsilon;
x<-seq(from=0,to=xMax,by=epsilon);

# generating the plot for the average of both groups
#title <- paste("B0 for Patient: ",patExp," from ", phenotype, " (blue)", "\n on ", arm, seg, " in 2D for ", dataSet, " data", sep="");
#title <- paste("Classifying patient (black): ",patExp);
#pdf(curvePath, width=9, height=6); # the initial height was 6
pdf(curvePath, width=9, height=3.5);
# par(mar= c(5.1, 4.1, 4.1, 2.1)) is the default
# bottom, left, top, right
par(mar= c(5.1, 4.1, 2.1, 2.1))
par(mfrow=c(1,2), cex.lab=1.2, cex.main=1.2);
#plot(x,AvgExp, type="p", pch=20, col="blue", ylim=c(0, yMax), xlab="Epsilon", ylab=param);
#points(x,AvgCon, pch=20, col="red");
plot(x,AvgExp, type="l" , lwd=3, col="blue", ylim=c(0, yMax), xlab="Epsilon", ylab=param);
lines(x,AvgCon, lty=5, lwd=3, col="red");
# now for the patient selected from test
points(x,data[patExp,], pch=20, col="black");
#title(title, adj=0);

# add a patient B0 curve
#title <- paste("B0 for Patient: ",patCon," from Non-", phenotype, " (red)","\n on ", arm, seg, " in 2D for ", dataSet, " data", sep="");
#title <- paste("Classifying patient (black): ",patCon);
plot(x,AvgExp, type="l", lwd=3, col="blue", ylim=c(0, yMax), xlab="Epsilon", ylab=param);
lines(x,AvgCon, lty=5,lwd=3, col="red");
# now for the patient selected from control
points(x,data[patCon,], pch=20, col="black");
#3title(title, adj=0);


dev.off();





