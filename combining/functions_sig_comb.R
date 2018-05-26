curvesMeansComb<-function(begPath, dataSet, param, dim, chrArm, phen1indices,phen2indices,seg=0){
  # Returns a data frame with two numeric vectors: test and control.
  # test and control contain the mean of conected components as epsilon increase
  # if we are running segments use seg=segment number
  # Note that this function is not working any more for parameter D
  
  
  # Read in data regardless of phenotype (jagged files from ~/Results/SET/2D/Homology)
  if (seg==0) {
    dataFile <- paste(param, "_", dim, "D_", dataSet, "_", chrArm, ".txt", sep="");
  } else {
    dataFile <- paste(param, "_", dim, "D_", dataSet, "_", chrArm,"_seg",seg, ".txt", sep="");
  }
  dataPath <- paste(begPath, "/Results/", dataSet, "/", subdir,"/",dim, "D/Homology/", dataFile, sep="");
  
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
