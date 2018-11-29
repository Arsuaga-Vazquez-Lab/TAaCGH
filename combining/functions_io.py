import csv, math, sys, os, random, errno

###################################
# Function to handle paths
# Thanks to Paul for the simplification

def makeDirectory(path):
	if os.path.isdir(path):
		return
	os.makedirs(path)

###################################
# Function to read from files

# A versatile reading function that casts numbers and can use
# a specified separator.
def readFile(rPath, separator, numCast):
	# Container for data to be returned
	rData = []

	rFile = open(rPath, "U")
	reader = csv.reader(rFile, dialect = 'excel', delimiter = separator)

	# Read the data
	for row in reader:
		tRow = []
		for j in range(len(row)):
			content = row[j]

			# Distinguish between strings versus type of number
			try:
				if numCast == "int":
					tRow.append(int(content))
				elif numCast == "float":
					tRow.append(float(content))
			except ValueError:
				tRow.append(content)
		rData.append(tRow)
	rFile.close()
	return rData

###################################
# Functions to write to files
###################################

# Write data to a separator delimited file
def writeFile (data, wPath, separator):
	# Get the directory from wPath
	wFilePath = os.path.dirname(wPath)

	try:
		# Check that the directory exists or can be created.
		makeDirectory(wFilePath)
		# Create the file object in write mode.
		dataFile = open(wPath,"w")
		
		#write to the file
		if type(data[0]) is list:
			for i in range(len(data)):
				for j in range(len(data[i])):
					dataFile.write(str(data[i][j]))
					if (j < len(data[i])-1):
						dataFile.write(separator)
				dataFile.write('\n')
			dataFile.close()
		else:
			for i in range(len(data)):
				dataFile.write(str(data[i]))
				dataFile.write('\n')
			dataFile.close()
	except OSError, ex:
		print >> sys.stderr, "There was an error in creating the path: %s" % (ex)
		raise ex

# For temporary files, we pass the fileObject and are careful not to close it.
def wTempTab (fileObject, data):
	#write to the file
	if type(data[0]).__name__ != 'list':
		for i in range(len(data)):
			fileObject.write(str(data[i]))
			fileObject.write('\n')
	else:
		for i in range(len(data)):
			for j in range(len(data[i])):
				fileObject.write(str(data[i][j]))
				if (j < len(data[i])-1):
					fileObject.write('\t')
			fileObject.write('\n')
	fileObject.seek(0)