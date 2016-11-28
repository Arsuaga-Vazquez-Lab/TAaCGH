# Each function which changes the data in some way will leave a receipt of the clones removed at each step

# This function filters a dataset based on the target information (could easily make it for the clones)
# We accept only targets/clones that are in the needles array. In needles array will come from a file
# which is a dictionary between clones, targets, chromosomes, and cytobands. For example, the UCSF HumArray 2.0
# file which can be downloaded from: http://0-www.ncbi.nlm.nih.gov.ilsprod.lib.neu.edu/projects/geo/query/acc.cgi?acc=GPL4420
filter.on.target <- function(cyto, needles, haystack, CGHstart) {
	print("Filtering on target...");
	nRowsNeedles <- nrow(needles);
		
	filteredHaystack <- c();
	
	for(i in c(1:nRowsNeedles)) {
		# Establish the needle, the indices in the haystack for that needle, and the data
		# in the haystack that goes with that needle.
		needleTarget <- toString(needles$Target[i]);
		haystackIndices <- which(as.matrix(haystack$Target) == needleTarget);
		haystackData <- haystack[haystackIndices, ];
		
		# We expect some Targets in haystack to not be in needles, but it is also the case that Targets
		# in needles could not be in haystack. This if accounts for that.
		if(length(haystackIndices) != 0) {
			# Now we're going to reconstitute the haystack with the same CGH information (obviously)
			# but with clone, chromosome, and base pair information from the UCSF HumArray 2.0 file
			# (needles) and arm/cytoband information from the UCSC Genome Browser (cyto).
			# For other datasets, we will still use the cyto file (possibly updated versions) and
			# we will use clone, target, chromosome, base pair information which describes whatever
			# array was used for the experiment.
			# This will standardize the data structure and order.
			
			# Grab what you can from the UCSF HumArray file.
			needleClone <- toString(needles$Clone[i]);
			needleChrom <- toString(needles$Chrom[i]);
			needleBP <- as.integer(needles$bp[i]);
			
			# Set up the grab from the cytoband file.
			# First, select the relevant chromosome information from the cytoband file.
			cytoIndices <- which(as.matrix(cyto$Chrom) == needleChrom);
			cytoData <- cyto[cytoIndices, ];
			nRowsCyto <- nrow(cytoData);
			
			# Determine extremes of the basepairs of the chromosome
			minBP <- as.integer(cytoData$Begin[1]);
			maxBP <- as.integer(cytoData$End[nRowsCyto]);
			
			# Make sure we'll find what we need. If the needleBP exceeds the maxBP from the
			# cytoBand file assign the arm and cytoband to the last row.
			if(needleBP > maxBP) {
				# stop(paste("The base pair of the needle (", needleClone,")", needleBP,"is outside [", minBP, ",", maxBP,"]", sep=" ", collapse=NULL));
				needleArm <- "q";
				needleCytoband <- toString(cytoData$Band[nRowsCyto]);
			} else {
				# Find the correct arm and cytoband
				for(j in c(1:nRowsCyto)) {
					# Check the row to see if needleBP is here and extract relevant information. 
					if(cytoData$Begin[j] < needleBP && needleBP < cytoData$End[j]) {
						needleArm <- toString(cytoData$Arm[j]);
						needleCytoband <- toString(cytoData$Band[j]);
						break;
					}
				}
			}
			
			# Stitch together the haystack data with the front matter.
			for(k in c(1:nrow(haystackData))) {
				row <- cbind(needleClone, needleTarget, needleChrom, needleArm, needleBP, needleCytoband, haystackData[k, CGHstart:ncol(haystackData)]);
				filteredHaystack <- rbind(filteredHaystack, row);
			}
		}
		
		# Print a progress report
		if(i %% 100 == 0) {
			print(paste("On row", i, "of", nRowsNeedles, sep=" ", collapse=NULL));
		}
	}
	
	# Name the front matter in a standard way. The remaining colnames are samples.
	colnames(filteredHaystack) <- c("Clone", "Target", "Chrom", "Arm", "bp", "Cytoband", colnames(haystack)[CGHstart:length(colnames(haystack))]);
	
	return(filteredHaystack);
}

# Assign the appropriate arm and cytoband based on the $bp column of needles, which here is the dataset
assign.arms <- function(cyto, needles, CGHstart) {
	print("Assigning arm information...");
	nRowsNeedles <- nrow(needles);
	
	armedNeedles <- c();
	
	for(i in c(1:nRowsNeedles)) {
		needleClone <- toString(needles$Clone[i]);
		needleChrom <- toString(needles$Chrom[i]);
		needleBP <- as.integer(needles$bp[i]);
		
		# Set up the grab from the cytoband file.
		# First, select the relevant chromosome information from the cytoband file.
		cytoIndices <- which(as.matrix(cyto$Chrom) == needleChrom);
		cytoData <- cyto[cytoIndices, ];
		nRowsCyto <- nrow(cytoData);
		
		# Determine extremes of the basepairs of the chromosome
		minBP <- as.integer(cytoData$Begin[1]);
		maxBP <- as.integer(cytoData$End[nRowsCyto]);
		
		# Make sure we'll find what we need. If the needleBP exceeds the maxBP from the
		# cytoBand file assign the arm and cytoband to the last row.
		if(needleBP > maxBP) {
			# stop(paste("The base pair of the needle (", needleClone,")", needleBP,"is outside [", minBP, ",", maxBP,"]", sep=" ", collapse=NULL));
			needleArm <- "q";
			needleCytoband <- toString(cytoData$Band[nRowsCyto]);
		} else {
			# Find the correct arm and cytoband
			for(j in c(1:nRowsCyto)) {
				# Check the row to see if needleBP is here and extract relevant information. 
				if(cytoData$Begin[j] <= needleBP && needleBP < cytoData$End[j]) {
					needleArm <- toString(cytoData$Arm[j]);
					needleCytoband <- toString(cytoData$Band[j]);
					break;
				}
			}
		}
		
		# Build the new data
		row <- cbind(needleClone, needleChrom, needleArm, needleBP, needleCytoband, needles[i, CGHstart:ncol(needles)]);
		armedNeedles <- rbind(armedNeedles, row);
		
		# Print a progress report
		if(i %% 100 == 0) {
			print(paste("On row", i, "of", nRowsNeedles, sep=" ", collapse=NULL));
		}
	}
	
	###################
	# This 2 really needs to be a 1 for the chin dataset because we get the bp information unlike with Climent...
	CGHstart <- CGHstart + 1;
	colnames(armedNeedles) <- c("Clone", "Chrom", "Arm", "bp", "Cytoband", colnames(armedNeedles)[CGHstart:length(colnames(armedNeedles))]);
	
	return(armedNeedles);
}

# Need to search for duplicate clones and average them to one clone
# Inputs are data (column names to be elaborated upon later), CGHstart (the column index at which CGH data starts)
# Will verify that clone names are the same and also check kb locations (throwing an error if they're not)
average.duplicates <- function(data, CGHstart) {
	print("Averaging duplicate probes...");
	
	nRows <- nrow(data);
	
	# Container for the result
	newData <- c();
	# Container for probes with duplicates and their number
	reportDropped <- c();
	
	i = 1;
	while(i < nRows) {
		# Define the needle, the haystack, and the indices where found
		cloneNeedle <- toString(data$Clone[i]);
		cloneHaystack <- as.matrix(data$Clone);
		foundIndices <- which(cloneHaystack == cloneNeedle);
		
		# Extract the duplicate clones
		duplicateClones <- data[foundIndices,];
		
		# If duplicateClones has more than one row, average over the columns while ignoring NAs.
		# Have to replace NaNs from entire column being NA with NA for later imputation.
		if(nrow(duplicateClones) != 1) {
			# Check the base pair locations are the same by seeing if standard deviation is 0
			genLocs <- as.numeric(as.matrix(duplicateClones$bp));
			if(sd(genLocs) != 0) {
				stop(paste("Error: Found clone match for", cloneNeedle,"but kb locations do not match.", sep=" ", collapse=NULL))
			}
			
			# Save the list of duplicate clones in their original form
			reportDropped <- rbind(reportDropped, duplicateClones);
			
			# Column means and replacement of NaN with NA
			averagedCGH <- colMeans(duplicateClones[,CGHstart:ncol(data)], na.rm=T);
			averagedCGH[is.nan(averagedCGH)] <- NA;
			
			# Put together the averaged data with the front matter information for the clone.
			# Grab the front matter of the first item in duplicateClones since it's the same for all.
			newData <- rbind(newData, c(duplicateClones[1,1:CGHstart-1], averagedCGH));
		} else {
			# If there are no duplicates just add the unadulterated entry to newData.
			newData <- rbind(newData, duplicateClones);
		}
		
		# Advance to the next row outside of the (possibly duplicate) block
		i = i + nrow(duplicateClones);
		
		# Print a progress report
		if(i %% 100 == 0) {
			print(paste("On row", i, "of", nRows, sep=" ", collapse=NULL));
		}
	}
	
	# Return both the new data containing the averages of duplicate clones and the
	# original clone data which were averaged.
	return(list(newData, reportDropped));
}

# Reduce according to proportion of missing NAs
row.reduce.NAs <- function(data, CGHstart, maxMissingProp) {
	print("Eliminating rows with too many NAs...");
	
	nRows <- nrow(data);
	
	reducedData <- c();
	reportDropped <- c();
	
	for(i in c(1:nRows)) {
		propNA <- length(which(is.na(data[i,CGHstart:ncol(data)]))) / (ncol(data) - CGHstart);
		if(propNA < maxMissingProp) {
			reducedData <- rbind(reducedData, data[i,]);
		} else {
			reportDropped <- rbind(reportDropped, data[i,]);
		}
		
		# Print a progress report
		if(i %% 100 == 0) {
			print(paste("On row", i, "of", nRows, sep=" ", collapse=NULL));
		}
	}
	
	return(list(reducedData, reportDropped));
}

# Verify the base pair orderings are correct are ascending within each chromosome
verify.kb.order <- function(data) {
	
}

# Verify that the clone ordering is the same
verify.clone.order <- function(data1, data2) {
	print("Verifying clone parity...");
	
	# Check that the lengths are the same
	if(nrow(data1) != nrow(data2)) {
		print("More clones in one dataset than the other. Apply match.clones.");
		return(list(FALSE, "Length"));
	}
	
	# At this point the lengths are the same so we should make
	# sure the order is actually the same.
	for(i in c(1:nrow(data1))) {
		if(toString(data1$Clone[i]) != toString(data2$Clone[i])) {
			# Pretty sure this never happens by virtue of the way match.clones works
			print("Clones are not in the same order. Apply order.clones.");
			return(list(FALSE, "Order"));
		}
	}
	
	# If we haven't left by now...
	print("The clones in the datasets are identical.");
	return(list(TRUE, "Great"));
}

# Match clones in two datasets. This will ensure that the set of clones in each dataset
# is the same, but will not necessarily guarantee the same ordering. That's for order.clones.
intersect.clones <- function(data1, data2) {	
	print("Intersect clones...");
	
	data1Clones <- as.matrix(data1$Clone);
	data2Clones <- as.matrix(data2$Clone);
	
	commonClones <- intersect(data1Clones, data2Clones);
	
	# The order here is CRUCIAL. When wrapping %in% in a which statement, the indices
	# returned are with respect to the X in X %in% Y. This is counter-intuitive, but a
	# simple example clears it up.
	data1Indices <- which(data1Clones %in% commonClones);
	data2Indices <- which(data2Clones %in% commonClones);
	
	# Extract the commonClones from each dataset
	matchedData1 <- data1[data1Indices, ];
	matchedData2 <- data2[data2Indices, ];
	
	# Report which clones have been removed from each dataset.
	droppedClones1 <- setdiff(data1Clones, commonClones);
	droppedClones2 <- setdiff(data2Clones, commonClones);
	
	droppedIndices1 <- which(data1Clones %in% droppedClones1);
	droppedIndices2 <- which(data2Clones %in% droppedClones2);
	
	# Extract the dropped clones from each dataset
	droppedData1 <- data1[droppedIndices1, ];
	droppedData2 <- data2[droppedIndices2, ];
	
	return(list(matchedData1, matchedData2, droppedData1, droppedData2));
}

# Order out-of-order data (oooData) based on a masterClones column vector. This will be called
# after matching is performed, so they will necessarily have the same length and the
# same membership!
order.clones <- function(masterClones, oooData) {
	print("Ordering clones to match the master list...");
	
	nRows <- nrow(masterClones);
	
	# The ordered data
	oData <- c();
	
	for(i in c(1:nrow(masterClones))) {
		# What are we trying to find?
		cloneNeedle <- toString(masterClones[i,1]);
		
		# Get the index of the clone in oooData.
		# We can be sure of a unique answer because all duplicates are gone.
		cloneIndexOoo <- which(oooData$Clone == cloneNeedle);
		
		oData <- rbind(oData, oooData[cloneIndexOoo, ]);
		
		# Print a progress report
		if(i %% 100 == 0) {
			print(paste("On row", i, "of", nRows, sep=" ", collapse=NULL));
		}
	}
	
	return(oData);
}
