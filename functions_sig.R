
#################################################
# FOR PERMUTATION TESTS

# Our typical sum of square of differences test statistic of the column means.
get_stat_diff_squared <- function(matrix1, matrix2)
{
	return(sum((colMeans(matrix1) - colMeans(matrix2))^2));
}

# The data on which this is performed is such that the rows are the patients
# and the columns are the epsilon increments (either of the networks or the 
# simplicial complexes) and the data are the parameter values. 
randomize <- function(df1, df2, fill, Replace_Subject=FALSE, Nitrs=10000)
{
	m1 <- as.matrix(df1);
	m2 <- as.matrix(df2);
	
	# Record dimension data
	n1 <- nrow(m1);
	n2 <- nrow(m2);
	N1 <- ncol(m1);
	N2 <- ncol(m2);
	
	# Fill in the missing rectangular matrix with the appropriate fill value so m1 and m2 are of the same size.
	# NOTE: In sig_pcalc.R, the data is already rectangularized so any block of patients has the same
	# number of columns as any other block of patients. This is here for generality.
	if(N1 > N2) {
		m2 <- cbind(m2,matrix(fill, nrow=n2, ncol=N1-N2));
		N2 <- N1;
	} else if(N2 > N1) {
		m1 <- cbind(m1,matrix(fill, nrow=n1, ncol=N2-N1));
		N1 <- N2;
	}
	
	if(N1 != N2) {
		stop("The number of columns does not match! The extension statements didn't work...");
	}
	
	# Compute the reference test statistic
	test_stat <- get_stat_diff_squared(m1,m2);
	# area by epsilon test-control (positive if blue on top)
	dif <- colMeans(m1) - colMeans(m2)
	# a column with a flag for positive and negative numbers (1, -1, 0)
	#sign<-sign(dif)
	# get the positive area
	#sumPos<-sum(dif[sign!=-1])
	sumPos<-sum(dif[sign(dif)==1])
	# get the negative area
	sumNeg<-sum(dif[sign(dif)==-1])
	
	if (test_stat==0) { 
		stop("The test statistic is 0! Something isn't right...");
	}
	
	# Combine the two classes of data
	temp_m <- rbind(m1,m2);
	
	# Create a vector of the row numbers
	x <- seq(1, nrow(temp_m));
	
	# Make a vector to store the test_stat for the permutations
	rnd_vec <- replicate(Nitrs,
							{
								# Create the permutated group indices
								group1 <- sample(x, n1, replace=Replace_Subject);
								group2 <- sample(x[-group1], n2, replace=Replace_Subject);
		
								# Calculate the statistic for the permutations
								get_stat_diff_squared(temp_m[group1,], temp_m[group2,]);
							}
							);

	return(list(sum(rnd_vec >= test_stat)/Nitrs, test_stat, sumPos, sumNeg));
}

########## FUNCTION curvesMeans

curvesMeans<-function(begPath, dataSet, subdir, param, dim, chrArm, phen1indices,phen2indices,seg=0){
  # Returns a data frame with two numeric vectors: test and control.
  # test and control contain the mean of conected components as epsilon increase
  # if we are running segments use seg=segment number
  # Note that this function is not working any more for parameter D
  
  
  # Read in data regardless of phenotype (jagged files from ~/Results/SET/2D/Homology)
  if (seg==0) {
    dataFile <- paste(param, "_", dim, "D_", dataSet, "_", subdir,"_",chrArm, ".txt", sep="");
  } else {
    dataFile <- paste(param, "_", dim, "D_", dataSet, "_",subdir,"_", chrArm,"_seg",seg, ".txt", sep="");
  }
  dataPath <- paste(begPath, "/Results/", dataSet, "/",subdir,"/", dim, "D/Homology/", dataFile, sep="");
  
  print(dataPath);
  
  # Need to determine the maximum number of columns in the file ahead of time
  # because, again, R is strange about reading in jagged arrays.
  nColData <- max(count.fields(dataPath, sep="\t"));
  
  # Read in parameter data
  data <- read.table(dataPath, header=FALSE, sep="\t", fill=TRUE, col.names=1:nColData);
  
  # Fill in the NAs in the jagged array with a parameter appropriate value.
  # For B0 and CC that value is 1 because that is the minimum/maximum, respectively.
  # For B1 that value is 0.

  if (identical(param,"B1")) {
    fillValue <- 0.0;
    data[is.na(data)] <- fillValue;
  } else if (identical(param,"CC") || identical(param,"B0")) {
    fillValue <- 1.0;
    data[is.na(data)] <- fillValue;
  } 
  
  phen1data <- data[phen1indices,];
  phen2data <- data[phen2indices,];
  
  phen1curve <- colMeans(phen1data);
  phen2curve <- colMeans(phen2data);
  
  curvesMeans=data.frame(test=phen1curve, control=phen2curve);	
  return(curvesMeans);
}

