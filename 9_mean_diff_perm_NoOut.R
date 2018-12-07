# Computing CENTERS OF MASS FOR EACH ARM and testing for difference (H0:mean(test)=mean(ctrl))
# This program computes permutations to get the p-value from the difference of means of a phenotype vs the rest for every full arm
# it also adjuts for multiple comparisons with FDR. 
# After averaging all clones in a chromosome arm, the center of mass is the average across the
# patients from either the test group or the control group. 
# Outliers registered at the phenotype file as "out_chrArm" (e.g. "out_17q") will
# be dropped before computing the center of mass of the arm

# INPUT
# 1) data set (SET_data_full.txt)
# 2) phenotype file (SET_phen.txt) with outliers by Chrom-Arm combination (if any)
# outliers should be registered by creating a column named "out_chrArm" (e.g. "out_17q")
# with 1 for an outlier and 0 for the rest of the patients. There is no need to have
# a column for a chrom-arm combination when there are not any outliers.

# OUTPUT
# Two files saved under ~/Reseach/Results/SET/CenterMass. They are  "SET_mean_diff_perm_phenotype_seedxx.txt" and
# "SET_mean_diff_perm_phenotype_seedxxsig.txt". The first file will contain the centers of mass for every arm
# The second file will be a subset with only those arms with FDR<sig (smaller FDR than the
# selected significance level). For the second file it will also determine if it is a gain or a loss
# or if it is undetermined using the decision tree below (flow chart)
# file containing: TestObsNum, ControlObsNum, mean_test,mean_ctrl, obs_mean_diff, sd_test, sd_ctrl, 
# t_test_0, t_ctrl_0, p_value_perm, fdr_perm
# In order of mention: Number of observations for test and control, mean for test group (phenotype chosen),
# mean for the rest of obervations (control), observed difference between mean_test and mean control, standard deviation for
# test group and control, t-test pvalue for center of mass (Ho: mean=0) for test phenotype and control phenotype, 
# p-value computed from permutations and fdr adjust for multiple comparisons

# ARGUMENTS
# 4. Name of dataset/file
# 5. segLength: Section size, use the same number used in the dictionary
#    This will be the minimum number of probes needed to run a specific arm
# 6. phenotype (lumA, rec, test, etc)
# permutations: number of permutations used to modify 0 p-values 
#                (e.g. 1/10,000=0.0001, pvalue[0]=0.0001)
# sig: desired false discovery rate (fdr) for the significant sections
# seed: an integer as a seed to have reproducible research

# EXAMPLE
# R --slave --args set 20 test 10000 0.05 1 < 9_mean_diff_perm_NoOut.R
# R --vanilla --args colon 50 younger45 10000 0.05 1 < 9_mean_diff_perm_NoOut.R

# Decision tree (flow chart)
# 1. Both probe averages (test & control) have the same sign or one of them is zero?
#     1.1 Yes (same sign). |test|>|ctrl|?
#         1.1.1 Yes. test>0?
#               1.1.1.1 test>0. it is a "GAIN" (Case 3.1)
#               1.1.1.2 test<=0 is is a "LOSS" (Case 3.2)
#         1.1.2 No. "Center of Mass driven by control aberrations" (Case 2.1)
#     1.2 No (different sign). Sample size for control >29 ?
#         1.2.1 Yes (nctrl>29). t-test H0:ctrl=0
#               1.2.1.1 Reject (H0:ctrl=0).  Sample size for test > 29?
#                        1.2.1.1. Yes (ntest>29). t-test H0: test=0
#                               1.2.1.1.1 Reject. "Undetermined: analyst must 
#                                         perform more tests to decide" (Case 7a.1)
#                               1.2.1.1.2 Not rejected.  "Center of Mass driven by 
#                                         control aberrations" (Case 7a.2)
#                        1.2.1.2 No (ntest<=29). "Undetermined: analyst must perform 
#                                                 more tests to decide" (Case 6a)
#               1.2.1.2 Not rejected (H0:ctrl=0). test>0?
#                       1.2.1.2.1 test>0. it is a "GAIN" (Case 5.1)
#                       1.2.1.2.2 test<0. it is a "LOSS" (Case 5.2)
#          1.2.2 No.(nctrl<=29). Sample size for test > 29?
#               1.2.2.1 Yes (ntest>29).  t-test H0: test=0
#                             1.2.2.1.1 Reject. "Undetermined: analyst must perform more 
#                                               tests to decide" (Case 7b.1)
#                             1.2.2.1.2 Not rejected.  "Center of Mass driven by 
#                                                   control aberrations" (case 7b.2)
#               1.2.2.1 No (ntest <=29). "Undetermined: analyst must perform more 
#                                         tests to decide" (Case 6b)

# Get the command line arguments
args = commandArgs();

dataSet <- args[4];
segLength<- as.numeric(args[5]);
phenotype <- args[6];
permutations <- as.numeric(args[7]);
sig <- as.numeric(args[8]);
seed <- as.numeric(args[9]);

# Only for debugging purposes
# dataSet <- "colon";
# segLength<-20;
# phenotype <- "younger45";
# permutations <- 10000;
# sig <- 0.05;
# seed <- 1;

# Working directory
#begPath <- "~/Research";
begPath <- "..";

# The following CGH_start is minus 1
CGH_start_minus1 <- 5;

### INPUT/OUTPUT FILES####

# Read the phenotype data
phenFile <- paste(dataSet, "phen.txt", sep="_");
phenPath <- paste(begPath, "Data", dataSet, phenFile, sep="/");
phenData <- read.table(phenPath, header=TRUE, sep="\t");

# Get the CGH data
dataPath <- paste(begPath, '/Data/', dataSet, '/', dataSet, '_data_full.txt', sep='');
data <- read.table(dataPath, header=T, sep='\t', comment.char='');

# Output directory
Folder = paste(begPath, 'Results', dataSet, 'CenterMass', sep='/');
if(!file.exists(Folder)) {
  dir.create(Folder, recursive=TRUE);
}


###############################
# BEGIN PROGRAM

# getting the indices for each phenotype group 
phen1indices <- which(phenData[,phenotype] == 1);
phen2indices <- which(phenData[,phenotype] == 0);
phen1num<-length(phen1indices);
phen2num<-length(phen2indices);

print(paste("phenotype", phenotype, sep=" "));

print(paste("There are ", phen1num, " patients for test", sep=""));
print(paste("There are ", phen2num, " patients for control", sep=""));

# container for results
res.t<-data.frame(phen1.num=integer(),phen2.num=integer(),m1=numeric(),m2=numeric(),mean_diff=numeric(),sd1=numeric(),sd2=numeric(),t_t=numeric(), t_c=numeric(), p.value=numeric());

# Determine the set of chromosomes in data set
chromList <- data$Chrom;
chroms <- unique(chromList);
arms <- c("p", "q");

# initializing variables
t.chr<-as.numeric();
t.arm<-as.character();

for(chr in chroms) {
  for(arm in arms) {

	print(paste("On chromosome ", chr, arm, sep=""));
	
	# Find the right arm
	dataIndices <- intersect(which(data$Chrom == chr), which(data$Arm == arm));
	profile_data <- data[dataIndices,];
	
	# keeping only arms with the minimum number of clones for the section
	if (length(dataIndices)>=segLength ) {
	  
	  # Creating a list of chromosom-arms
	  t.chr<-rbind(t.chr,chr);
	  t.arm<-rbind(t.arm,arm);

    # Remove outliers from arm
	  # Initializing variables
    phen1indices_noOut<-phen1indices;
    phen2indices_noOut<-phen2indices;
    phen1num_noOut<-phen1num;
    phen2num_noOut<-phen2num;
    outlier<-paste("out_",chr,arm, sep="");
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
      
    # Get patients average profile for phenotype
    # Working on phen1
    phen1_profile_data <- profile_data[, (CGH_start_minus1 + phen1indices_noOut)];
	  print("phen1");
	  print(phen1_profile_data[1,c(1:5)]);
	  ave_phen1 <- colMeans(phen1_profile_data);
	  # Descriptive statistics 
	  gran_ave_phen1 <- mean(ave_phen1);
	  sd_phen1 <- sd(ave_phen1);
      
    # Working on phen2
	  phen2_profile_data <- profile_data[, (CGH_start_minus1 + phen2indices_noOut)];
	  print("phen2");
	  print(phen2_profile_data[1,c(1:5)]);
	  ave_phen2 <- colMeans(phen2_profile_data);
	  # Descriptive statistics
	  gran_ave_phen2 <- mean(ave_phen2);
	  sd_phen2 <- sd(ave_phen2);
	
	  obs_mean_diff<-gran_ave_phen1-gran_ave_phen2;
	  
	  # t-test, one sample, mu=0, two-sided
	  t_test_0<-t.test(ave_phen1, mu=0);
	  t_ctrl_0<-t.test(ave_phen2, mu=0);
	
	  # merging both phenotypes for permutation 
	  merged<-cbind(phen1_profile_data,phen2_profile_data);
	
	  # Permutation
	  # Create a vector of the row numbers
	  index <- seq(1, phen1num_noOut+phen2num_noOut);
	
	  # Using a seed for reproducible research
	  set.seed(seed);	
	
	  # Make a vector to store the test_stat for the permutations
	  perm_mean_diff <- replicate(permutations, {
		  										# Create the permutated group index just for the first group 
			  									# (the rest is the other group)
				  								sample <- sample(index, phen1num_noOut, replace=FALSE);
					 
					  							# Calculate the statistic for the permutations: difference of means
						  						mean(colMeans(merged[,sample]))-mean(colMeans(merged[,-sample]));
							  				}
								  );
	  # Computing the p-value from the permutations
	  p_value_perm<-sum(abs(perm_mean_diff)>=abs(obs_mean_diff))/permutations;
	  res.t<-rbind(res.t, c(phen1num_noOut, phen2num_noOut, gran_ave_phen1, gran_ave_phen2, obs_mean_diff, sd_phen1, sd_phen2, t_test_0$p.value, t_ctrl_0$p.value, p_value_perm));
 	  }
   }
 }

colnames(res.t)<-c("TestObsNum","ControlObsNum","mean_test","mean_ctrl","obs_mean_diff","sd_test", "sd_ctrl", "t_test_0", "t_ctrl_0", "p_value_perm");

# the following option is to avoid scientific notation in the output
options("scipen"=100);

# Avoiding 0 p-value before adjusting with FDR
noZero<-1/permutations;
res.t$p_value_perm[res.t$p_value_perm==0]<-noZero

# Adjust for multiple testing with FDR
res.t$fdr_perm<-round(p.adjust(res.t$p_value_perm,"fdr"), digits=4);

#Adding chromosome and arm information
res.t$Chr<-t.chr;
res.t$Arm<-t.arm;

# Formatting output
res.t$mean_test<-round(res.t$mean_test, digits=2);
res.t$mean_ctrl<-round(res.t$mean_ctrl, digits=2);
res.t$obs_mean_diff<-round(res.t$obs_mean_diff, digits=2);
res.t$sd_test<-round(res.t$sd_test, digits=2);
res.t$sd_ctrl<-round(res.t$sd_ctrl, digits=2);
res.t$t_test_0<-round(res.t$t_test_0, digits=4);
res.t$t_ctrl_0<-round(res.t$t_ctrl_0, digits=4);
res.t$p_value_perm<-round(res.t$p_value_perm, digits=4);

# output file name and path
File <- paste(dataSet, '_mean_diff_perm_noOut_', phenotype, '_seed',seed, sep='');
Path <- paste(Folder, "/", File, ".txt", sep='');

write.table(res.t, file=Path, row.names=F, sep='\t')

######## SELECTING AND SAVING SIGNIFICANT SECTIONS ###########
##### DETERMINING IF CM FOR ARMS IS GAIN OR LOSS ####
sig_fdr <- subset(res.t, res.t$fdr<=sig);

# Initializing Arm_conclusion
sig_fdr$Arm_conclusion<-rep(NA, nrow(sig_fdr));

if (nrow(sig_fdr)>0) {
  for (i in 1:nrow(sig_fdr)) {
    print("i="); print(i); 
    # 1. If mean_test and mean_ctrl have the same sign
    if (sig_fdr$mean_test[i]*sig_fdr$mean_ctrl[i]>=0) {
      print("1. same sign");
      # 2. AND |mean_test|>|mean_ctrl| then keep significance
      if (abs(sig_fdr$mean_test[i])>abs(sig_fdr$mean_ctrl[i])) {
        print("2.|mean_test|>|mean_ctrl|?yes");
        # 3. Positive CM or Negative CM
          if (sig_fdr$mean_test[i]>0) {print("3.mean_test>0?yes"); 
                                      sig_fdr$Arm_conclusion[i]<-"Gain (Case 3.1)"
          } else {print("3.mean_test>0?no"); 
            sig_fdr$Arm_conclusion[i]<-"Loss (Case 3.2)";
          }# end of else from (3)
        
      # else from (2): |mean_test|<|mean_ctrl|
      } else {print("2.|mean_test|>|mean_ctrl|?no");
              sig_fdr$Arm_conclusion[i]<-"Center of Mass driven by control aberrations (Case 2.1)";
      } # end of else from (2)
    
    # else from (1): mean_test and mean_ctrl have the different sign
    } else {
      print("Different Sign");
      # 4. if sample size from control > 29 H0:mean_ctrl=0 t-test large sample
      if (sig_fdr$ControlObsNum[i]>29) {print("4.nctrl>29?yes");
        # 5. if H0 rejected: mean_ctrl might be different than 0
        if (sig_fdr$t_ctrl_0[i]<sig) {print("5.t_ctrl_0 rejected?yes");
          # 6a. if sample size from test > 29 H0:mean_test=0 t-test large sample
          if (sig_fdr$TestObsNum[i]>29) {
            # 7a. if H0 rejected: mean_test might be different than 0
            if (sig_fdr$t_test_0[i]<sig) {
              sig_fdr$Arm_conclusion[i]<-"Undetermined: analyst must perform more tests to decide (Case 7a.1)";
            # else from (7a): H0 not rejected, no evidence that mean_test<>0 (=0)
            } else sig_fdr$Arm_conclusion[i]<-"Center of Mass driven by control aberrations (Case 7a.2)";
          # else from (6a): sample size from test <= 29
          } else {
            sig_fdr$Arm_conclusion[i]<-"Undetermined: analyst must perform more tests to decide (Case 6a)";
          } # end of else from (6a)
        # else from (5): H0 not rejected, no evidence that mean_ctrl<>0 (=0)
        } else {print("5.t_ctrl_0 rejected?no (=0)");
          if (sig_fdr$mean_test[i]>0) {print("mean_test>0"); 
                                      sig_fdr$Arm_conclusion[i]<-"Gain (Case 5.1)"; 
          } else {print("mean_test<=0"); 
            sig_fdr$Arm_conclusion[i]<-"Loss (Case 5.2)";
          } # end of else from above
        } # end of else from (5)
      # else from (4): sample size from control <= 29
      } else {print("4.nctrl>29?no");
         # 6b. if sample size from test > 29 H0:mean_test=0 t-test large sample
         if (sig_fdr$TestObsNum[i]>29) {print("6b.ntest>29?yes");
            # 7b. if H0 rejected: mean_test might be different than 0
            if (sig_fdr$t_test_0[i]<sig) {print("7b.t_test_0 rejected?yes");
             sig_fdr$Arm_conclusion[i]<-"Undetermined: analyst must perform more tests to decide (Case 7b.1)";
            # else from (7b): H0 not rejected, no evidence that mean_test<>0 (=0)
            } else {print("7b.t_test_0 rejected?no (=0)");
              sig_fdr$Arm_conclusion[i]<-"Center of Mass driven by control aberrations (Case 7b.2)";
            } # end of else from (7b)
         # else from (6b): sample size from test <= 29
         } else {print("6b.ntest>29?no")
           sig_fdr$Arm_conclusion[i]<-"Undetermined: analyst must perform more tests to decide (Case 6b)";
         } # end of else from (6b)
      } # end of else from (4)
    } # end of else from (1)
  
  } # end of "for loop"  
      
  # Saving into file significant FDR sections
  sig_fdrName <- paste(File,"sig.txt", sep="");
  sig_fdrPath <- paste(Folder, sig_fdrName, sep="/");
  
  print(sig_fdrPath);
  write.table(sig_fdr,sig_fdrPath, sep='\t', row.names=FALSE);
# else from if (nrow(sig_fdr)>0)
} else print("There are not significant arms for displacement of center of mass")



