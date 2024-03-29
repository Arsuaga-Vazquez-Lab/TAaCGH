#################################################
# USAGE
#
# Given a dataset, this code will compute the betti numbers for the point clouds created
# from the aCGH data from either all full arms or all sections of every arm and by patient.
# To speed up the calculations 2_cgh_dictionary_cytoband.R splits the the work into parts
# that can be run simultaneously.
#
# We utilize code from the Comptop group at Stanford (jplex) to compute the homology:
# http://comptop.stanford.edu/u/programs/jplex/index.html.
#
# Parameters: The epsilon increments for the filtration, epsIncr. The part number which
# refers to the proper dictionary file, partNum.
#
# Input: 1) The tab delimited transposed version of aCGH data arranged as mentioned in readme.docx.
# The dictionary file (from cgh_dictionary.R) which indicates the 0-based
# start and stop indices of the chromosome arms 2) Dictionary files 
#
# Output:
# It creates two different kind of files:
# 1) Barcode files named Inter_2D_hom1_SET_ChromArm_patientNo_segNo.txt
# under ~/Research/Results/dataSet/action/2D/Homology/Chromosome_number  (action: arms or sect)
# 2) Jagged files (B0_2D_SET_ChromArm_segNo.txt) with the value for the homology for each patient (row)
# at each increment for epsilon (columns). The Epsilon increment is saved for future reference in a
# text file named Epsilon_XXX.txt. This file and all the jagged files will be saved under:
# ~/Research/Results/SET/2D/Homology/B0_2D_SET_chrArm_seg.txt

import csv, math, sys, os, tempfile, subprocess, re, string
from functions_io import *
from functions_cgh import *

#################################################
# COMMAND LINE ARGUMENTS
# 1. dataSet (sj, climent, sim, etc.)
# 2. homDim (usually 1 or 2, which computes B0 or (B0 and B1), respectively.
# 3. partNum (a positive integer)
# 4. epsIncr (usually 0.01 or 0.05)
# 5. action (arms or sect) specify if you are running full arms or sections (as used in 2_cgh_dictionary_cytoband.R)
# note: dataSet will be expected within Research/dataSet but diccionaries within a subfolder named arms or sect

# 
# EXAMPLE
# > python 4_hom_stats_parts.py set 1 2 0.01 action

dataSet = sys.argv[1]
homDim = int(sys.argv[2])
partNum = int(sys.argv[3])
epsIncr = round(float(sys.argv[4]), 3)
action = sys.argv[5]

#################################################
# FILE SPECIFIC FUNCTIONS

def makeDirectory(path):
    if os.path.isdir(path):
        return
    os.makedirs(path)

#################################################
# BEGIN PROGRAM

# Determine base path
#begPath = os.path.join(os.getenv('HOME'), "Research")
# will take any starting point instead of HOME
current = os.getcwd()
begPath=os.path.abspath(os.path.join(current, os.pardir))

# Set the paths for the dictionary, and data files
dictFile = "%s_%s_dict_%d.txt" % (dataSet, action, partNum)
#dictPath = os.path.join(begPath, "Data",dataSet, dictFile)
dictPath = os.path.join(begPath, "Data", dataSet, action, dictFile)

dataFile = "%s_data.txt" % (dataSet)
dataPath = os.path.join(begPath, "Data", dataSet, dataFile)

# Read data and dictionary
inputList = readFile(dataPath, "\t", "float")
dictList = readFile(dictPath, "\t", "int")

# Get rid of the first two rows.
inputList.pop(0)
inputList.pop(0)
dictList.pop(0)

# Using dimension 2 as window
cloudDim = 2
 
# Establish the results directory
resultsPath = os.path.join(begPath, "Results", dataSet, action, repr(cloudDim)+"D", "Homology")

# Go through each row of the dictionary list to get all chromosome and arm combinations.
for row in dictList:
	
	# Chromosome and arm assignments
    chr = row[0]
    arm = row[1]
	
	# Beginning and ending assignments
    beg = row[2]
    end = row[3]
	
	# Segment number
    seg = row[5]
	
    # Store the Betti numbers for all patients in a list containing as many lists as homDim - 1.
    bettiNums = [[] for d in range(homDim)]
		
    # Go through the individual samples.
    for m in range(len(inputList)):
		
        print "(Homology) On chromosome %d arm %s segment %d dimension %d and patient %d of %d" % (chr,arm, seg, cloudDim, (m+1), len(inputList))
            
        #################################################
        # Generate data for use in jPlex and run it.
        #################################################
			
        # Build the individual profile.
        profile = build_individual_profile(inputList, m, beg, end)
			
        # Build the cloud.
        cloud = build_cloud(profile, cloudDim)
			
        # Build the distance matrix from the cloud.
        distanceMatrix = build_distance_matrix(cloud, cloudDim)
			
        # Select as the largest epsilon, the largest distance divided by 3.
        maxEps = max(map(max, distanceMatrix)) * 0.9
			
        # Create the temporary file object which we'll reference now and again.
        cloudFile = tempfile.NamedTemporaryFile(suffix = '.txt')
        try:
            # Write the cloud to the temporary file.
            wTempTab(cloudFile, cloud)
            
            # File to write the intervals. Writing is handled by beanshell.
            makeDirectory(resultsPath + "/" + repr(chr))
            intFile = "Inter_%dD_hom%d_%s_%d%s_pat%d_seg%d.txt" % (cloudDim, homDim, dataSet, chr, arm, (m+1), seg)
            intPath = os.path.join(resultsPath, repr(chr), intFile)
            
            # When you have time, try going back and using tempfile for this. Weird that you couldn't get it to work the first time.
            # Create the bash file. The names will be unique so we can run in parallel with no problem.
            bashName = "bash%d%s%d%d%d%d%s.bsh" % (chr, arm, beg, cloudDim, m, seg, dataSet)
            bashFile = open(bashName, 'w')
				
            bashFile.write('pcloud = Plex.EuclideanArrayData("'+cloudFile.name+'");\n')
            bashFile.write('rips = Plex.RipsStream('+repr(epsIncr)+', '+repr(homDim)+', '+repr(maxEps)+', pcloud);\n')
            bashFile.write('intervals = Plex.Persistence().computeIntervals(rips);\n')
            bashFile.write('file = new FileOutputStream("'+intPath+'");\n')
            bashFile.write('p = new PrintStream(file);\n')
            bashFile.write('this.interpreter.setOut(p);\n')
            bashFile.write('print(intervals);')
				
            bashFile.close()
            
            # Run jPlex with the temporary cloud and bash files. Need to wait() or we may delete the bash file before it runs.
            plex = subprocess.Popen("java -Xmx1024m -cp plex.jar JTerm < "+bashFile.name, shell=True)
            plex.wait()
				
            # Remove the temporary bash file when we're done with it.
            os.remove(bashName)
        finally:
            # Automatically clean up the cloud file
            cloudFile.close()
			
        #################################################
        # Parse the interval output from jPlex
        #################################################
			
        # Open the homFile
        homFile = open(intPath, 'r')
			
        # Create regular expression object.
        # ((?<=\[)\d) gives the dimension of the interval.
        # (\d\.\d+) gives interval bounds (or starting point if it's inf)
        pattern = re.compile("((?<=\[)\d)|(\d\.\d+)")
			
        # Going to hold the cleaned up intervals
        intervals = []
			
        # Keep track of the maximum filtration time present
        maxFilt = 0.0
        
        # Go through each line of homFile and create an array: Dimension, Interval Start, Interval End
        for line in homFile:
            # Get the results of findall
            result = pattern.findall(line)
				
            tRow = []
				
            # If the result is non-empty, we're going to want to get the stuff.
            if len(result) != 0:
                # This will always be the dimension
                tRow.append(int(result[0][0]))
					
                # This will always be the start of the interval
                tRow.append(float(result[1][1]))
					
                # This is the end of the interval
                if len(result) == 3:
                    tRow.append(float(result[2][1]))
						
                    # Update maximum filtration time
                    if float(result[2][1]) > maxFilt:
                        maxFilt = float(result[2][1])
                else:
                    tRow.append('inf')
					
                intervals.append(tRow)
			
        # Close homFile because we're done reading it. Delete file located at intPath (messy version).
        homFile.close()
        os.remove(intPath)
			
        # Write to file located at intPath (clean version). 
        writeFile(intervals, intPath, "\t")
			
        #################################################
        # Collect incremental betti's from intervals
        #################################################
			
        # Store the betti numbers. This creates a list with homDim - 1 empty lists contained to hold the betti numbers.
        indBettiNums = [[] for d in range(homDim)]
			
        # Increment epsilons at the same incrementation when you made the bash script.
        eps = epsIncr
			
        # Go through eps = 0.005 until one increment after maxFilt
        while (eps <= maxFilt + epsIncr):
            # Store the temporary count of betti numbers
            tBettis = [0 for d in range(homDim)]
				
            # Go through the entries of the interval table and count if the epsilon is in the interval.
            for entry in intervals:
                if repr(entry[2]) == 'inf':
                    tBettis[entry[0]] = tBettis[entry[0]] + 1
                elif entry[1] <= eps and eps <= entry[2]:
                    tBettis[entry[0]] = tBettis[entry[0]] + 1
				
            # Append the counts for this epsilon
            for i in range(len(indBettiNums)):
                indBettiNums[i].append(tBettis[i])
				
            eps = eps + epsIncr
			
        # Append the patient's list of betti counts to the list containing all the patients (bettiNums)
        for i in range(len(indBettiNums)):
            bettiNums[i].append(indBettiNums[i])
		
    # Write the betti counts for all the patients in this chromosome arm to a file.
    for hDim in range(len(bettiNums)):
        bettiPath = "%s/B%d_%dD_%s_%s_%d%s_seg%d.txt" % (resultsPath, hDim, cloudDim, dataSet, action, chr, arm, seg)
        writeFile(bettiNums[hDim], bettiPath, "\t")

# Saving in a text file the epsilon increment used for the homology
EpsPath = "%s/Epsilon_%04.2f.txt" % (resultsPath, epsIncr)
EpsText = "Used epsilon=%04.2f" % (epsIncr)
makeDirectory(resultsPath)
EpsFile = open(EpsPath,"w")
EpsFile.write(EpsText)
EpsFile.close()
