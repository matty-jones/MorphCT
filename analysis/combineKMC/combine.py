import cPickle as pickle
import os
import subprocess as sp
import sys
from scipy.sparse import lil_matrix



if __name__ == "__main__":
    outputDir = sys.argv[1]
    tasksetID = sys.argv[2]
    currentPID = os.getpid()
    affinityJob = sp.Popen(['taskset', '-pc', str(tasksetID), str(currentPID)], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    print "Combining outputs..."
    combinedData = {}
    pickleFiles = []
    for fileName in os.listdir(outputDir):
        if ("KMCResults_" in fileName):
            pickleFiles.append(outputDir + '/' + fileName)
    pickleFiles = sorted(pickleFiles)
    print "%d pickle files found to combine!" % (len(pickleFiles))
    for fileName in pickleFiles:
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
    print "Writing out the combined pickle file..."
    with open(outputDir + '/KMCResults.pickle', 'w+') as pickleFile:
        pickle.dump(combinedData, pickleFile)
    print "Complete data written to", outputDir + "/KMCResults.pickle."
    print "Cleaning up..."
    # Delete any unneeded files
    sp.Popen("rm -f " + outputDir + "/KMCData*", shell = True, stdout = open(os.devnull, 'wb'), stderr = open(os.devnull, 'wb')).communicate()
    #sp.Popen("rm -f " + outputDir + "/KMCResults_*", shell = True, stdout = open(os.devnull, 'wb'), stderr = open(os.devnull, 'wb')).communicate()
