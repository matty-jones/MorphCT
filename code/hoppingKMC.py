import os
import sys
import numpy as np
import random as R
from scipy.sparse import lil_matrix
import cPickle as pickle
import subprocess as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import helperFunctions
try:
    import mpl_toolkits.mplot3d.axes3d as p3
except ImportError:
    print "Could not import 3D plotting engine, calling the plotCarrier subroutine will result in an error!"


def execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList):
    # Attempt 2. PROGRAM THIS SERIALLY FIRST SO THAT IT WORKS
    # Determine the maximum simulation times based on the parameter dictionary
    simulationTimes = parameterDict['simulationTimes']
    carrierList = []
    # Create the carrierList which contains the information that the singleCore program needs to run its jobs
    for carrierNo in range(parameterDict['numberOfCarriersPerSimulationTime']):
        for lifetime in simulationTimes:
            carrierList.append([carrierNo, lifetime])
    # The carrierList is now like the ORCAJobsList, so split it over each procID
    procIDs = parameterDict['procIDs']
    outputDir = parameterDict['outputDir'] + '/' + parameterDict['morphology'][:-4] + '/KMC'
    jobsList = [carrierList[i:i + (int(np.ceil(len(carrierList) / len(procIDs)))) + 1] for i in xrange(0, len(carrierList), int(np.ceil(len(carrierList)/float(len(procIDs)))))]
    print "Writing job pickles for each CPU..."
    runningJobs = []
    for procID, jobs in enumerate(jobsList):
        pickleName = outputDir + '/KMCData_%02d.pickle' % (procID)
        with open(pickleName, 'w+') as pickleFile:
            pickle.dump(jobs, pickleFile)
        print "KMC jobs for procID", procID, "written to KMCData_%02d.pickle" % (procID)
        # Open the required processes to execute the KMC jobs
        print 'python ' + os.getcwd() + '/code/singleCoreRunKMC.py ' + outputDir + ' ' + str(procID) + ' &'
        runningJobs.append(sp.Popen(['python', str(os.getcwd()) + '/code/singleCoreRunKMC.py', outputDir, str(procID)]))
    # Wait for all jobs to complete
    [p.wait() for p in runningJobs]
    # Now combine all of the pickle files into one:
    print "All KMC jobs completed!"
    if parameterDict['combineKMCResults'] is True:
        print "Combining outputs..."
        combinedData = {}
        for procID, jobs in enumerate(jobsList):
            fileName = outputDir + '/KMCResults_%02d.pickle' % (procID)
            # The pickle was repeatedly dumped to, in order to save time.
            # Each dump stream is self-contained, so iteratively unpickle to add the new data.
            with open(fileName, 'r') as pickleFile:
                pickledData = pickle.load(pickleFile)
                for key, val in pickledData.iteritems():
                    if key not in combinedData:
                        combinedData[key] = val
                    else:
                        combinedData[key] += val
        # Write out the combined data
        with open(outputDir + '/KMCResults.pickle', 'w+') as pickleFile:
            pickle.dump(combinedData, pickleFile)
        print "Complete data written to", outputDir + "/KMCResults.pickle."
        print "Cleaning up..."
        # Delete any unneeded files
        sp.Popen("rm -f " + outputDir + "/KMCResults_*", shell = True, stdout = open(os.devnull, 'wb'), stderr = open(os.devnull, 'wb')).communicate()
    sp.Popen("rm -f " + outputDir + "/KMCData*", shell = True, stdout = open(os.devnull, 'wb'), stderr = open(os.devnull, 'wb')).communicate()
    return AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList

if __name__ == "__main__":
    try:
        pickleFile = sys.argv[1]
    except:
        print "Please specify the pickle file to load to continue the pipeline from this point."
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = helperFunctions.loadPickle(pickleFile)
    execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList)
