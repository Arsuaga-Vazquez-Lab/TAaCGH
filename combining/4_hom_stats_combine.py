
# coding: utf-8

# In[2]:

import csv, math, sys, os, tempfile, subprocess, re, string
from functions_io import *
from functions_cgh import *


# In[ ]:

#################################################
# COMMAND LINE ARGUMENTS
# 1. dataSet (sj, climent, sim, etc.)
# 2. homDim (usually 1 or 2, which computes B0 or (B0 and B1), respectively.
# 3. partNum (a positive integer)
# 4. epsIncr (usually 0.01 or 0.05)
# 5. subdir the name of a subdir within /dataSet dir for this particular run
# note: dataSet will be expected within /dataSet dir but diccionaries within subdir
# results for homology will be saved under /Results/dataSet

# EXAMPLE
# > python 4_hom_stats_parts.py set 1 2 0.01 subdir

dataSet = sys.argv[1]
homDim = int(sys.argv[2])
partNum = int(sys.argv[3])
epsIncr = round(float(sys.argv[4]), 3)
subdir = sys.argv[5]


# In[3]:

# In case of running from the console
#dataSet = 'kwek8p11q'
#homDim = 1
#partNum = 1
#epsIncr = 0.1
#subdir = 'arms'


# In[4]:


#################################################
# FILE SPECIFIC FUNCTIONS

def makeDirectory(path):
    if os.path.isdir(path):
        return
    os.makedirs(path)


# In[5]:

#################################################
# BEGIN PROGRAM

# Determine base path
begPath = os.path.join(os.getenv('HOME'), "Research")

# Reading corresponding dictionary
dictFile = "%s_dict_comb_%d.txt" % (dataSet, partNum) 
dictPath = os.path.join(begPath, "Data", dataSet, subdir, dictFile)
dictList = readFile(dictPath, "\t", "int")
# Get rid of header
dictList.pop(0)


# In[6]:

# Reading data
dataFile = "%s_data.txt" % (dataSet)
dataPath = os.path.join(begPath, "Data", dataSet, dataFile)
inputList = readFile(dataPath, "\t", "float")
# Get rid of header (2 lines)
inputList.pop(0)
inputList.pop(0)


# In[20]:

# Using dimension 2 as window
cloudDim = 2

# Establish the results directory
resultsPath = os.path.join(begPath, "Results", dataSet, subdir, repr(cloudDim)+"D", "Homology")
resultsPath


# In[22]:

# Go through each row of the dictionary list to get all chromosome and arm combinations.
for row in dictList:
	
	# info from first section (arms)
    chr1 = row[0]
    arm1 = row[1]
    beg1 = row[2]
    end1 = row[3]
	
	# info from second section (arms)
    chr2 = row[4]
    arm2 = row[5]
    beg2 = row[6]
    end2 = row[7]
	
	
    # Store the Betti numbers for all patients in a list containing as many lists as homDim - 1.
    bettiNums = [[] for d in range(homDim)]
		
    # Go through the individual samples.
    for m in range(len(inputList)):
		
        print "(Homology) On combination %d %s %d %s and patient %d of %d" % (chr1,arm1,chr2,arm2,(m+1), len(inputList))
            
        #################################################
        # Generate data for use in jPlex and run it.
        #################################################
			
        # Build the individual profile from first section.
        profile = build_individual_profile(inputList, m, beg1, end1)
			
        # Build the cloud first section.
        cloud = build_cloud(profile, cloudDim)
        
        # Build the individual profile from second section.
        profile = build_individual_profile(inputList, m, beg2, end2)
        
        # Build the cloud second section.
        cloud2 = build_cloud(profile, cloudDim)
        
        # combine both clouds
        #cloud = cloud.append(cloud2)
        cloud += cloud2
			
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
            makeDirectory(resultsPath + "/" + repr(chr1))
            intFile = "Inter_%dD_hom%d_%s_%d%s_%d%s_pat%d.txt" % (cloudDim, homDim, dataSet, chr1, arm1,chr2, arm2,(m+1))
            intPath = os.path.join(resultsPath, repr(chr1), intFile)
            
            # When you have time, try going back and using tempfile for this. Weird that you couldn't get it to work the first time.
            # Create the bash file. The names will be unique so we can run in parallel with no problem.
            bashName = "bash%d%s%d_%d%s%d_%d%d%s.bsh" % (chr1, arm1, beg1,chr2, arm2, beg2, cloudDim, m, dataSet)
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
        bettiPath = "%s/B%d_%dD_%s_%d%s_%d%s.txt" % (resultsPath, hDim, cloudDim, dataSet, chr1, arm1,chr2, arm2)
        writeFile(bettiNums[hDim], bettiPath, "\t")



# In[23]:

# Saving in a text file the epsilon increment used for the homology
EpsPath = "%s/Epsilon_%04.2f.txt" % (resultsPath, epsIncr)
EpsText = "Used epsilon=%04.2f" % (epsIncr)
makeDirectory(resultsPath)
EpsFile = open(EpsPath,"w")
EpsFile.write(EpsText)
EpsFile.close()


# In[ ]:



