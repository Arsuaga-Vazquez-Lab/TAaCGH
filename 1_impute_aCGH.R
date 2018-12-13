# TODO this script was used to impute some of our datasets but i hasn't been updated since 2016
# This program impute missing values from aCGH data with lowess using the
# aCGH package from Bioconductor
# (http://www.bioconductor.org/packages//2.11/bioc/html/aCGH.html)

# INPUT
# aCGH file. It must be a tab delimited text file (.txt)
# the file must follow the format as specified in README.txt
# Use a short name for the file, lets say "set" and save it under ~/Research/Data/set

# OTHER REQUIREMENTS
# the library aCGH from Bioconductor must be installed in advance

# Arguments
# 1. fileName without including the txt extension (e.g. "set")

# OUTPUT
# a file named set_lowess.txt will be saved in ~/Research/Data/set

# Example of the command in the terminal (use vanilla instead of
# slave to have input from the terminal)
# R --slave --args set < 1_impute_aCGH.R

# Get the command line arguments
args = commandArgs();
fileName<-args[4];


library(aCGH);


begPath <- "~/Research";
CGH_start <- 6


# Read the data
dataPath <- paste(begPath, "/Data/", fileName, "/", fileName, ".txt", sep = "");
data <- read.table(dataPath, sep = "\t", header = TRUE);

# Format log2.ratio file from data so aCGH package can be used
log2ratios <- data[, CGH_start: ncol(data)]; 
rownames(log2ratios) <- data$Clone;
 
# create clones.info for the aCGH package: mapping information including but not
# limited to clone name, chromosome and kb relative to the chromosome
clones_info <- data[, c(1:(CGH_start -1))];
# Change bp to kb
clones_info$bp <- round(clones_info$bp / 1000);
names(clones_info)[names(clones_info)=="bp"]<-"kb"

# Create aCGH object
ex_acgh <- create.aCGH(log2ratios, clones_info);

# Impute missing data
log2ratios_imputed <- impute.lowess(ex_acgh);
# put imputed data back together with the front matter...
lowess_imputed_data <- cbind(data[,1:(CGH_start -1)], log2ratios_imputed);

# Write lowess imputed climent data to file
write.table(lowess_imputed_data, paste(begPath, "/Data/", fileName, "/",fileName, "_lowess.txt", sep=''), sep='\t', row.names=F);