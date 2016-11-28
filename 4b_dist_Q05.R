# This program computes for each segment the minimum distance between points
# aswell as 5% percentile and then take the average across patients
# Avg_Min and Avg_Q05 in the dictionary

# R --vanilla --args dataSet < dist_Q05.R
# R --vanilla --args bergamaschi1p < dist_Q05.R

# Get the command line arguments
args = commandArgs();

dataSet <- args[4];

# Only for debugging
# dataSet<-"horlings_sect"

begPath <- "~/Research";
dim<-"2";

###############################
# READ FILES

# Read the complete chromosome dictionary data
chrDictFile <- paste(dataSet, "_dict_cyto.txt", sep="");
chrDictPath <- paste(begPath, "Data", dataSet, chrDictFile, sep="/");
chrDict <- read.table(chrDictPath, header=TRUE, sep="\t");
print(chrDictPath);

# Get the CGH data
dataPath <- paste(begPath, '/Data/', dataSet, '/', dataSet, '_data_full.txt', sep='');
data <- read.table(dataPath, header=T, sep='\t', comment.char='');


#################
# BEGIN PROGRAM

Avg_Min<-c()
Avg_Q05<-c()
# Reading dictionary to go segment by segment
for (i in 1:nrow(chrDict)) {
  chr <- chrDict[i,1];
  arm <- as.vector(chrDict[i,2]);
  beg<-chrDict[i,3];
  end<-chrDict[i,4];
  chrDict[i,4];
  seg <- chrDict[i,6];
  print(paste("On chromosome ", chr,arm, " segment ", seg, sep=""));
  
  # Initializing
  minAllPats<-c();
  sum<-0;
  
  # cloud and minimum euclidean distance
  for (pat in 6:ncol(data)) {
    cloud<-cbind(data[(beg+1):(end+1),pat], c(data[(beg+2):(end+1),pat], data[(beg+1),pat]))
    dist<-dist(cloud, method="euclidean")
    minAllPats<-c(minAllPats,min(dist))
    sum<-sum+quantile(dist,0.05)
  }
  Avg_Min<-c(Avg_Min,mean(minAllPats))
  Avg_Q05<-c(Avg_Q05,sum/(ncol(data)-5))
}
chrDict$Avg_Min<-Avg_Min
chrDict$Avg_Q05<-Avg_Q05

# Let's save these statistics in the dictionary   
write.table(chrDict,chrDictPath, sep='\t', row.names=FALSE)
    