import csv, math, sys, os

# Grab the CGH profile for a particular arm.
def build_individual_profile(CGHdata , patient , start , finish):
	tProfile = CGHdata[patient][start:(finish + 1)]
	return tProfile

# Build the n-dimensional cloud of points from a CGH profile.
def build_cloud(CGHprofile , dim):
	tempPoint = []
	tCloud = []
	# Go through the number of points for that chromosome arm
	for i in range(len(CGHprofile)):
		# Have as many coordinates as dimensions
		for j in range(dim):
			# This ensures proper wraparound.
			tempPoint.append(round(CGHprofile[((i+j)%(len(CGHprofile)))],6))
		tCloud.append(tempPoint)
		tempPoint = []
	return tCloud

# Build the distance matrix from the sorted cloud. ith row jth column.
# Append 0.0001 as the distance if points are identical so that we'll
# always consider them connected from the start, and we'll avoid infinite
# loops for duplicate points. This is arguably much simpler than finding
# duplicate points and deleting them.
def build_distance_matrix(pointCloud , dim):
	tSumOfSquaresOfDiffs = 0
	tDistMat = []
	tRow = []
	for i in range(len(pointCloud)):
		for j in range(len(pointCloud)):
			if i < j:
				tSumOfSquaresOfDiffs = 0
				for k in range(dim):
					tSumOfSquaresOfDiffs += (pointCloud[j][k]-pointCloud[i][k])**2
				tDist = math.sqrt(tSumOfSquaresOfDiffs)
				if tDist != 0:
					tRow.append(tDist)
				else:
					tRow.append(0.0001)
			else:
				tRow.append(0)
		tDistMat.append(tRow)
		tRow = []
	return tDistMat
	
# Convert a distance matrix to an adjacency matrix given a distance threshhold
def build_adjacency_matrix(distMat, dist):
	adjMat = []
	adjRow = []
	
	for row in distMat:
		for elem in row:
			if elem < dist and elem != 0:
				adjRow.append(1)
			else:
				adjRow.append(0)
		adjMat.append(adjRow)
		adjRow = []
	
	return adjMat