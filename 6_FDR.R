# TODO do I need to add I way to run a specific directory without running all the parts (like when you need to run a subset of the dictionary)? Need a parameter to tell the program to run a specific dictionary
# This program will adjust the p-values from 5_sig_pcalc_parts.R for multiple comparisons using FDR
# It will read the different files (parts) and will create a single file with adjusted p-value
# This script will need to be run only once no matter how many parts were used to split the dictionary
# "Idx0.Beg" and "Idx0.End" are the beginning and ending positions in the CGH file for the sections (0 start)

# INPUT
# p-value files, if the data set was called SET they must be in ~Research/Results/SET/significance/pvals
# pvalue file names look like B0_sim_SET_pvals_1.txt (sim=phephenotype, SET=dataset and 1=part number)
# Notice that the Parts Argument refers to the TOTAL number of parts

# OUTPUT
# 2 files: one ending in xxx_FDR.txt with all the sections and their corresponding False Dsicovery Rate (FDR)
# the other file ending in xxx_FDRsig.txt with a subset of sections from xxx_FDR.txt with FDR<=sig
# where sig is the threshold chosen as the maximum FDR allowed for significance.

# ARGUMENTS
# 4. Name of dataset/file
# 5. Parameter (B0, B1, etc)
# 6. phenotype (lumA, test, sim, rec, etc)
# 7. Parts (in which the directory was split (TOTAL))
# 8. perm: number of permutations used to modify 0 p-values (e.g. 1/10,000=0.0001, pvalue[0]=0.0001)
# 9. sig: desired false discovery rate (fdr) for the significant sections

# EXAMPLE
# R --slave --args set B0 test 2 10000 0.05< 6_FDR.R
# R --vanilla --args bergamaschiMADMA3_sect B0 Luminal_A 7 10000 0.05< 6_FDR.R

# Get the command line arguments
args = commandArgs();
 
file <- args[4];
param <- args[5];
phen <- args[6];
parts <- args[7];
perm <- as.numeric(args[8]);
sig <- as.numeric(args[9]);

# Only for debugging purposes
# file <- "bergamaschi1pMADMA3";
# param <- "B0";
# phen <- "ErbB2";
# parts <- 1;
# perm <- 10000;
# sig <- 0.05;
 
# Working directory
begPath <- "~/Research";

######## CREATING A SINGLE P-VALUE FILE FROM THE PARTS ###########
# Read the first file
begName <- paste(param, phen, file, "pvals", sep="_");
first <- paste(begName, "1.txt", sep="_");
Path <- paste(begPath, "Results", file, "significance", "pvals", sep="/");
filePath <- paste(Path, first, sep="/");
print(filePath);
all_fdr <- read.table(filePath, header=TRUE, sep="\t");
all_fdr<-all_fdr[order(all_fdr$Chr, all_fdr$Arm,all_fdr$Segment), ];

# Paste together all the files
if (parts>1) {
  for (i in c(2:parts)) {
	  tempPath <- paste(Path, "/", begName, "_", i, ".txt", sep=""); 
	  print(tempPath);
	
	  tempFile <- read.table(tempPath, header=TRUE, sep="\t");
	  tempFile<-tempFile[order(tempFile$Chr, tempFile$Arm,tempFile$Segment), ];
	  all_fdr <- rbind(all_fdr, tempFile);
	  rm(tempFile);
  }
}

######## ADJUSTING FOR MULTIPLE TESTING WITH FDR ########

# Avoiding 0 p-value before adjusting with FDR
noZero<-1/perm;

all_fdr$P.Value<-ifelse(all_fdr$P.Value==0, noZero, all_fdr$P.Value);

# Adjust for multiple testing with FDR
all_fdr$fdr<-round(p.adjust(all_fdr$P.Value,"fdr"), digits=4);

# Formating output for Test.Stat and P.Value
all_fdr$Test.Stat<-round(all_fdr$Test.Stat, digits=2);
all_fdr$P.Value<-round(all_fdr$P.Value, digits=4);

#Saving the complete fdr adjusted p-values into file
all_fdrName <- paste(begName, "FDR", sep="_");
all_fdrPath <- paste(Path, "/", all_fdrName,".txt", sep="");

print(all_fdrPath);

write.table(all_fdr,all_fdrPath, sep='\t', row.names=FALSE);

######## SELECTING AND SAVING SIGNIFICANT SECTIONS ###########
sig_fdr <- subset(all_fdr, all_fdr$fdr<=sig);

# Saving into file significant FDR sections
sig_fdrName <- paste(all_fdrName,"sig.txt", sep="");
sig_fdrPath <- paste(Path, sig_fdrName, sep="/");

print(sig_fdrPath);
write.table(sig_fdr,sig_fdrPath, sep='\t', row.names=FALSE)

