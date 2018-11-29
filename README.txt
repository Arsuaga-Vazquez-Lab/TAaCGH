Readme

This text describes the steps to use TAaCGH software. Each main script starts with a number signaling the order to run them.
Take a look into the code of each script as it provides useful information about the arguments to
use as well as the input and output. There are also other pieces of code for which their names
do not start with a number. Those are either functions called by the scripts or jplex.

1. Create a directory called ~/Research/Code and save all the programs there. All the programs are
meant to be run in the terminal by calling a command with the necessary arguments. Information about
the arguments and examples of the command lines are available in the first lines of the code.
2. Pre-process aCGH file. Lets call the tab delimited file "SET.txt"
	2.1 Update the position (bp) for every clone using the last build from a genome browser.
	2.2 The file should have 5 columns with the information about the clone and its position. 
	The names of those columns must be "Clone" "Chrom" "Arm" "Cytoband" "bp". 
	For X chromosome use 23 and for Y use 24 in humans (otherwise will produce weird ordering to dataSet). 
	aCGH information for every patient(row) shall start in column number 6.
	2.3 Verify that there are no duplicates in the file for every Chrom/Arm/bp 
	combination. If there are any, they must be averaged
	2.4 Impute missing values using impute_aCGH.R. This program uses lowess from the
	aCGH package from Bioconductor (http://www.bioconductor.org/packages//2.11/bioc/html/aCGH.html)
	2.5 Name: Use a short name for the file, for instance "SET" 
	and add at the end "data_full.txt" so that the name of the file become something like 
	SET_data_full.txt and save it under ~/Research/Data/SET	
3. Phenotype file. Should be named SET_phen.txt. Consists of two or more columns "patID" and one or more columns for the different phenotypes. Each phenotype column must be coded with 1 for those patient in the test group and 0 for those in control. patID will have the ID of the patient (do not start the ID with a number, make sure the ID follows the rules for variable names in R). The program works using an index, therefore it is very important that the file is ordered by patient in the exact same order as the patients in the aCGH file. And of course it should have the exact same patients, no more, no less.
3.1 Inspect the data for outliers and find the chromosome and arm for those clones. If any, register them by chrom-arm combination creating a column in the phenotype file named "out_chrArm" (e.g. "out_17q") with 1’s for outliers and 0’s for the rest of the patients. There is no need for a column if the chrom-arm combination did not have any outliers. The patient will not be considered when computing the center of mass of that particular arm
4. From SET_data_full.txt, create a "dictionary" specifying the index positions for each segment within the arms 
using 2_cgh_dictionary_cytoband.R. Decide the number of parts to split the file to speed up the the homology and the permutations by running each part in parallel. This script will overwrite set_data_full.txt after reordering the 
dataSet and will generate set_data_orig.txt with the original dataSet.
4b. Compute summary statistics (Average minimum and average 5% percentile) for the distance between points using 3b_dist_Q05.R. Results will be saved in SET_dict_cyto.txt
5. Create SET_data.txt, the transposed version of SET_data_full.txt using 3_Transposed_aCGH.R
6. Run the homology (B0 or B1) feeding SET_data.txt into 4_hom_stats_parts.py. You will need to run it as many times as the number of parts used to split the dictionary
7. Compute the un-adjusted p-values using 5_sig_pcalc_parts.R. You will need to run it as many times as the number of parts used to split the dictionary
8. Adjust the p-values using 6_FDR.R. You will need to run this script only once even though it uses all parts with the un-adjusted p-values
9. Generate B0 curves for significant sections in xxx_FDRsig.txt (from 6_FDR.R) using 7_vis_curves.R. Inspect manually the output and keep only as significant those sections with the test curve (blue) above the control curve (red). Generate a reviewed file with
only the new significant sections and name it xxxx_FDRsig_rev.txt
10. Find the significant probes from significant sections using 8_probesFDR.R feeding xxx_FDRsig_rev.txt into it. 
The output is two files, one with p-values from all significant sections and a second file with only those probes under the desired
significant level. This second file will also show if the probe is a gain or a loss
11. 9_mean_diff.perm.R This program computes the centers of mass for each chromosome-arm combination and performs the statistical test for the difference between test and control for a specific phenotype. The  only input is the aCGH file SET_data_full.txt and the phenotype file. It could be run after step 3. As in step 10 there are two output files, one with all centers of mass and a second one displaying
only those arms with a p-value under the significance value and it also shows if it is a gain or a loss

Requirements for the software:
aCGH library from R: used in 1_impute_aCGH.R (http://www.bioconductor.org/packages//2.11/bioc/html/aCGH.html)
Python to run  4_hom_stats_parts.py
JDX (Java) needed to run jPlex in 4_hom_stats_parts.py

Improvements:
+ Handling outliers for Centers of Mass: 05/05/2016
