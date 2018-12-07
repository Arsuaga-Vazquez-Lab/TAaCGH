# This program uses permutations to find significant probes (Ho: test=control) within a
# significant segment from 6_FDR.R and 7_vis_curves.R

# INPUT
# a manually reviewed file from xxxx_FDRsig.txt under the name xxxx_FDRsig_rev.txt. The file should contain
# those significant sections that have the B0 curve for test above the control curve (7_vis_curves.R)

# OUTPUT
# two files will be saved under ~Results/SET/significance/pvals, they are xxx_Probes_FDR.txt
# and xxx_Probes_FDRsig.txt. The first file will contain all the p-values and FDRs from the input
# file and the second will be a subset with only those clones with FDR<sig (smaller FDR than the
# selected significance level). For the second file it will also determine if it is a gain or a loss
# or if it is undetermined using the decision tree from below (flow chart)

# ARGUMENTS
# 4. Parameter (B0, B1, etc)
# 5. phenotype (lumA, test, sim, rec, etc)
# 6. Name of dataset/file
# 7. subdir (where the dictionaries were saved and where the p-values were saved)
# 8. perm: number of permutations used to modify 0 p-values (e.g. 1/10,000=0.0001, pvalue[0]=0.0001)
# 9. sig: desired false discovery rate (fdr) for the significant sections
# 10. seed: any integer number for reproducible research

# EXAMPLE
# R --slave --args B0 test set subdir 10000 0.05 1 < 8_probesFDR.R
# R --vanilla --args B0 TP53_mut horlings sect 10000 0.05 1 < 8_probesFDR.R

# Decision tree (flow chart)
# 1. Both probe averages (test & control) have the same sign  or one of them is zero?
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
#                                                 perform more tests to decide" (Case 7a.1)
#                               1.2.1.1.2 Not rejected.  "Center of Mass driven by 
#                                                         control aberrations" (Case 7a.2)
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
#                                                       control aberrations" (case 7b.2)
#               1.2.2.1 No (ntest <=29). "Undetermined: analyst must perform more 
#                                         tests to decide" (case 6b)


# Get the command line arguments
args = commandArgs();

param <- args[4];
phen <- args[5];
file <- args[6];
subdir <- args[7];
perm <- as.numeric(args[8]);
sig <- as.numeric(args[9]);
seed <- as.numeric(args[10]);

# Only for debugging purposes
# file <- "horlings_sect";
# param <- "B0";
# phen <- "ER_pos";
# perm <- 10000;
# sig <- 0.05;
# seed <- 1;

# Working directory
#begPath <- "~/Research";
begPath <- "..";

# CGH_start
CGH_start_minus1 <- 5;

####### Input files ########

# Read the file with significant sections
begName <- paste(param, phen, file, subdir, "pvals_FDRsig_rev", sep="_");
first <- paste(begName, ".txt", sep="");
resPath <- paste(begPath, "Results", file, subdir,"significance", "pvals",param,phen, sep="/");
sections_Path <- paste(resPath, first, sep="/");
sig_sections <- read.table(sections_Path, header=T, sep='\t', comment.char='');

# Read the phenotype data
phenFile <- paste(file, "phen.txt", sep="_");
phenPath <- paste(begPath, "Data", file, phenFile, sep="/");
phenData <- read.table(phenPath, header=TRUE, sep="\t");

# Get the CGH data
dataPath <- paste(begPath, '/Data/', file, '/', file, '_data_full.txt', sep='');
CGH <- read.table(dataPath, header=T, sep='\t', comment.char='');

###### Output file name and path #######
probes_fdr <- paste(param, "_", phen, "_", file, '_Probes_FDR.txt', sep='');
Path <- paste(resPath, probes_fdr, sep='/');

# the following option is to avoid scientific notation in the output
options("scipen"=100);

###############################
# BEGIN PROGRAM

# getting the indices for each phenotype group 
phen1indices <- which(phenData[,phen] == 1);
phen2indices <- which(phenData[,phen] == 0);
phen1num<-length(phen1indices);
phen2num<-length(phen2indices);

print(paste("phenotype", phen, sep=" "));
print("phenData$phenotype", phenData[,phen]);

print(paste("There are ", phen1num, " patients for test", sep=""));
print(paste("There are ", phen2num, " patients for control", sep=""));

# subset CGH data to probes within significant sections from the FDRsig file
# initializing
beg <- sig_sections$Idx0.Beg[1]+1;
end <- sig_sections$Idx0.End[1]+1;
subCGH<-CGH[c(beg:end),];

if (nrow(sig_sections)>1) {
  # Cycle through the FDRsig file
  for(i in c(1:nrow(sig_sections))) {
    # Get chromosome and arm
    beg=ifelse(sig_sections$Idx0.Beg[i]+1<end,end+1,sig_sections$Idx0.Beg[i]+1)
    end <- sig_sections$Idx0.End[i]+1;
    subCGH<-rbind(subCGH,CGH[c(beg:end),]);
  }
}
####### The real stuff for the test #####

# Get average probe by phenotype
phen1_data<-subCGH[, (CGH_start_minus1 + phen1indices)]; 
print("phen1");
print(phen1_data[1,c(1:5)]);
avg_probe_phen1<-round(rowMeans(phen1_data), digits=2);
phen2_data<-subCGH[, (CGH_start_minus1 + phen2indices)]; 
print("phen2");
print(phen2_data[1,c(1:5)]);
avg_probe_phen2<-round(rowMeans(phen2_data), digits=2);

# Descriptive statistics
sd_probe_phen1<-round(apply(subCGH[, (CGH_start_minus1 + phen1indices)],1,sd), digits=2);
sd_probe_phen2<-round(apply(subCGH[, (CGH_start_minus1 + phen2indices)],1,sd), digits=2);
obs_mean_diff<-avg_probe_phen1-avg_probe_phen2;

# t-test, one sample, mu=0, two-sided
t_test_0<-numeric(nrow(subCGH));
t_ctrl_0<-numeric(nrow(subCGH));
for (i in 1:nrow(subCGH)) {
  t_test<-t.test(phen1_data[i, ], mu=0);
  t_test_0[i]<-round(t_test$p.value, digits=4);
  t_ctrl<-t.test(phen2_data[i, ], mu=0);
  t_ctrl_0[i]<-round(t_ctrl$p.value, digits=4);
}

# merging both phenotypes for permutation 
merged<-cbind(phen1_data,phen2_data);

# Permutation
# Create a vector of the col numbers
index <- seq(1, phen1num+phen2num);

# Using a seed for reproducible research
set.seed(seed);	

# Permutations: Matrix storing in each column the test_stat {avg(probes.test)-avg(probes.ctrl)} for each permutation (rows=probes)
perm_mean_diff <- replicate(perm, 
                            {
                              # Create the permutated group index just for the first group 
                              # (the rest is the other group)
                              sample <- sample(index, phen1num, replace=FALSE);
                              
                              # compute the average per row (probe) by phenotype and substract them (code by probe)
                              rowMeans(merged[,sample])-rowMeans(merged[,-sample]);
                            }
);

p_value_perm<-rowSums(abs(perm_mean_diff)>=abs(obs_mean_diff))/perm;
p_value_perm[p_value_perm==0]<-1/perm;
perm_fdr<-round(p.adjust(p_value_perm,"fdr"),digits=4);

# output
subCGH <- data.frame(lapply(subCGH[,1:CGH_start_minus1], as.character), stringsAsFactors=FALSE);
rows<-length(subCGH$Clone);
testnum<-rep(phen1num,rows);
ctrlnum<-rep(phen2num,rows);
output<-cbind(testnum,ctrlnum,subCGH,avg_probe_phen1,avg_probe_phen2,round(obs_mean_diff, digits=3),sd_probe_phen1,sd_probe_phen2,t_test_0, t_ctrl_0, round(p_value_perm, digits=4),perm_fdr);
colnames(output)<-c("TestObsNum","ControlObsNum","Clone","Chrom","Arm","bp","Cytoband","avg_probe_test","avg_probe_ctrl","obs_mean_diff","sd_probe_test","sd_probe_ctrl","t_test_0", "t_ctrl_0", "p_value_perm","perm_fdr");

write.table(output, file=Path, row.names=F, sep='\t')

######## SELECTING AND SAVING SIGNIFICANT PROBES ###########
# Saving into file significant FDR sections
sig_fdr <- subset(output, output$perm_fdr<=sig);

# Initializing Arm_conclusion
sig_fdr$Arm_conclusion<-rep(NA, nrow(sig_fdr));

if (nrow(sig_fdr)>0) {
  for (i in 1:nrow(sig_fdr)) {
    print("i="); print(i); 
    # 1. If avg_probe_test and avg_probe_ctrl have the same sign or one is zero
    if (sig_fdr$avg_probe_test[i]*sig_fdr$avg_probe_ctrl[i]>=0) {
      print("1. same sign");
      # 2. AND |avg_probe_test|>|avg_probe_ctrl| then keep significance
      if (abs(sig_fdr$avg_probe_test[i])>abs(sig_fdr$avg_probe_ctrl[i])) {
        print("2.|avg_probe_test|>|avg_probe_ctrl|?yes");
        # 3. Positive CM or Negative CM
        if (sig_fdr$avg_probe_test[i]>0) {print("3.avg_probe_test>0?yes"); 
                                          sig_fdr$Arm_conclusion[i]<-"Gain (Case 3.1)"
        } else {print("3.avg_probe_test>0?no"); 
                sig_fdr$Arm_conclusion[i]<-"Loss (Case 3.2)";
        } # end of else from (3)
      # else from (2): |avg_probe_test|<|avg_probe_ctrl|
      } else {print("2.|avg_probe_test|<|avg_probe_ctrl|?no");
            sig_fdr$Arm_conclusion[i]<-"Center of Mass driven by control aberrations (Case 2.1)";
      } # end of else from (2)
    
    # else from (1): avg_probe_test and avg_probe_ctrl have the different sign
    } else { print("Different Sign");
      # 4. if sample size from control > 29 H0:avg_probe_ctrl=0 t-test large sample
      if (phen2num>29) {print("4.nctrl>29?yes");
        # 5. if H0 rejected: avg_probe_ctrl might be different than 0
        if (sig_fdr$t_ctrl_0[i]<sig) {print("5.t_ctrl_0 rejected?yes");
          # 6a. if sample size from test > 29 H0:avg_probe_test=0 t-test large sample
          if (phen1num>29) {
            # 7a. if H0 rejected: avg_probe_test might be different than 0
            if (sig_fdr$t_test_0[i]<sig) {
              sig_fdr$Arm_conclusion[i]<-"Undetermined: analyst must perform more tests to decide (Case 7a.1)";
              # else from (7a): H0 not rejected, no evidence that avg_probe_test<>0 (=0)
            } else sig_fdr$Arm_conclusion[i]<-"Center of Mass driven by control aberrations (Case 7a.2)";
            # else from (6a): sample size from test <= 29
          } else {
            sig_fdr$Arm_conclusion[i]<-"Undetermined: analyst must perform more tests to decide (Case 6a)";
          } # end of else from (6a) 
        # else from (5): H0 not rejected, no evidence that avg_probe_ctrl<>0 (=0)
        } else {print("5.t_ctrl_0 rejected?no (=0)"); 
          if (sig_fdr$avg_probe_test[i]>0) {print("avg_probe_test>0");
                                      sig_fdr$Arm_conclusion[i]<-"Gain (Case 5.1)";
          } else {print("avg_probe_test<=0");
                  sig_fdr$Arm_conclusion[i]<-"Loss (Case 5.2)";
          } # end of else from above
        } # end of else from (5)
      # else from (4): sample size from control <= 29
      } else {print("4.nctrl>29?no");
        # 6b. if sample size from test > 29 H0:avg_probe_test=0 t-test large sample
        if (phen1num>29) {print("6b.ntest>29?yes");
          # 7b. if H0 rejected: avg_probe_test might be different than 0
          if (sig_fdr$t_test_0[i]<sig) {print("7b.t_test_0 rejected?yes");
            sig_fdr$Arm_conclusion[i]<-"Undetermined: analyst must perform more tests to decide (Case 7b.1)";
          # else from (7b): H0 not rejected, no evidence that avg_probe_test<>0 (=0)
          } else {print("7b.t_test_0 rejected?no (=0)");
            sig_fdr$Arm_conclusion[i]<-"Center of Mass driven by control aberrations (Case 7b.2)";           } # end of else from (7b)
        # else from (6b): sample size from test <= 29
        } else {print("6b.ntest>29?no");
          sig_fdr$Arm_conclusion[i]<-"Undetermined: analyst must perform more tests to decide (Case 6b)";
        }  # end of else from (6b)
      } # end of else from (4)
    } # end of else from (1)
  
  } # end of "for loop"  
  
  sigProbes_fdr <- paste(param, "_", phen, "_", file, '_Probes_FDRsig.txt', sep='');
  sig_fdrPath <- paste(resPath, sigProbes_fdr, sep='/');
  print(sig_fdrPath);
  write.table(sig_fdr,sig_fdrPath, sep='\t', row.names=FALSE);
# else from if(nrow(sig_fdr>0))
} else print("There are not significant clones within the significant sections")
