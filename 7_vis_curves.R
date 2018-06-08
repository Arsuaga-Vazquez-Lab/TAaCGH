# This program generate B0 curves for the significant sections resultant
# from 6_FDR listed under xxxx_FDRsig.txt and will save them under 
# ~Research/Results/dataSet/vis/curves/2D/B0/phenotype

# INPUT
# 1. The phenotype file where the patient names matches the order from the aCGH file with the log2 ratios.
# A column with the phenotype (ERBB2, basal, test, sim, etc) coded with 0's and 1's. The column name will be the
# argument "phenotype"
# 2. The jagged array (from 4_hom_stats_parts.py) containing the value of the parameter (B0, B1) for all the patient
# simplicial complexes/networks over the epsilon filtrations. The files should be under 
# ~/Results/SET/2D/Homology/B0_2D_SET_chrArm_seg.txt for a specific aCGH file named SET and a specific chromosome-arm
# combination "chrArm"


#################################################
# COMMAND LINE ARGUMENTS
# 4. param (B0, B1)
# 5. phenotype (ERBB2, basal, test, sim, etc)
# 6. dataSet (SET, etc)
# 7. subdir (where the dictionaries where saved)

# EXAMPLE
# R --slave --args B0 test set < 7_vis_curves.R
# R --slave --args B0 TP53_mut horlings_sect < 7_vis_curves.R

# Get the command line arguments
args = commandArgs();

param <- args[4];
phenotype <- args[5];
dataSet <- args[6];
subdir <- args[7];

# for debugging purposes only
# param <- "B1";
# phenotype <- "both8p11q";
# dataSet <- "kwek8p11q";
# subdir <-  "sects";
   
###############################
# READ FILES

# Establish the beginning path
begPath <- "~/Research";

srcPath <- paste(begPath, "Code", "functions_sig.R", sep="/");
source(srcPath);

# Read the file with significant sections
begName <- paste(param, phenotype, dataSet, subdir,"pvals_FDRsig", sep="_");
# use the following line to see all curves (not only significant)
#begName <- paste(param, phenotype, dataSet, subdir,"pvals","FDR", sep="_");
Path <- paste(begPath, "Results", dataSet, subdir, "significance", "pvals", sep="/");
filePath <- paste(Path, "/", begName,".txt", sep="");
print(filePath);
sig_sections <- read.table(filePath, header=TRUE, sep="\t");

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

# Ensure the folder we are going to write to exists
curveFolder = paste(begPath, 'Results', dataSet, subdir, 'vis', 'curves', "2D", param, phenotype, sep='/');
if(!file.exists(curveFolder)) {
dir.create(curveFolder, recursive=TRUE);
}


# Go through the rows of the file with significant sections
for(i in c(1:nrow(sig_sections)))
{
	chr <- sig_sections$Chr[i];
	arm <- as.vector(sig_sections$Arm[i]);
	chrArm <- paste(chr, arm, sep="");
	seg <- sig_sections$Segment[i];
	
	print(paste("On chromosome ", chrArm, " segment ", seg,  sep=""));
	
	# Plots will be saved under the following name
	curveFile <- paste(param, '_', phenotype, '_2D_', dataSet, '_', chrArm, '_s', seg, '.pdf', sep='');
	curvePath <- paste(curveFolder, curveFile, sep='/');
	
	
	# generating the ith curveMeans
	curvesMeans_i<-curvesMeans(begPath=begPath, dataSet=dataSet, param=param, dim=2, chrArm=chrArm, phen1indices=phen1indices, phen2indices=phen2indices, seg=seg);
	phen1curve <- curvesMeans_i$test
	phen2curve <- curvesMeans_i$control
	yMax <- max(max(phen1curve), max(phen2curve));
	xMax<-length(phen1curve)*0.01-0.01;
	x<-seq(from=0,to=xMax,by=0.01);
	
	# generating the ith plot
	title <- paste(param, " curves for ", phenotype, " (", phen1num, " blue) vs Non-", phenotype, " (", phen2num, " red) on ", chr, arm, ".s", seg, " in 2D for ", dataSet, " data", sep="");
	pdf(curvePath, width=9, height=6);
	par(mfrow=c(1,1), cex.lab=1.2, cex.main=1.2);
	plot(x,phen1curve, type="p", pch=20, col="blue", ylim=c(0, yMax), xlab="Epsilon", ylab=param);
	points(x,phen2curve, pch=20, col="red");
	title(title, adj=0);
	dev.off();
}
	
	
	

