import os
import glob
import sys
import numpy as np
import random as R
from scipy.sparse import lil_matrix
import pickle
import subprocess as sp
from morphct.definitions import SINGLE_RUN_MOBKMC_FILE
from morphct.code import helperFunctions


def execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList):
    try:
        if parameterDict['useAverageHopRates']:
            print("Be advised: useAverageHopRates is set to", repr(parameterDict['useAverageHopRates']) + ".")
            print("ORCA-calculated energy levels will be ignored, and the following hop rates will be used:")
            print("Average Intra-molecular hop rate:", parameterDict['averageIntraHopRate'])
            print("Average Inter-molecular hop rate:", parameterDict['averageInterHopRate'])
    except KeyError:
        pass
    # Attempt 2. PROGRAM THIS SERIALLY FIRST SO THAT IT WORKS
    # Determine the maximum simulation times based on the parameter dictionary
    simulationTimes = parameterDict['simulationTimes']
    carrierList = []
    # Modification: Rather than being clever here with the carriers, I'm just going to create the master
    # list of jobs that need running and then randomly shuffle it.
    # This will hopefully permit a similar number of holes and electrons and lifetimes to be run simultaneously
    # providing adequate statistics more quickly
    for lifetime in simulationTimes:
        for carrierNo in range(parameterDict['numberOfHolesPerSimulationTime']):
            carrierList.append([carrierNo, lifetime, 'Hole'])
        for carrierNo in range(parameterDict['numberOfElectronsPerSimulationTime']):
            carrierList.append([carrierNo, lifetime, 'Electron'])
    R.shuffle(carrierList)
    # Old method:
    ## Create the carrierList which contains the information that the singleCore program needs to run its jobs
    #for carrierNo in range(parameterDict['numberOfHolesPerSimulationTime']):
    #    for lifetime in simulationTimes:
    #        carrierList.append([carrierNo, lifetime, 'Hole'])
    #for carrierNo in range(parameterDict['numberOfElectronsPerSimulationTime']):
    #    for lifetime in simulationTimes:
    #        carrierList.append([carrierNo, lifetime, 'Electron'])
    # The carrierList is now like the ORCAJobsList, so split it over each procID
    procIDs = parameterDict['procIDs']
    outputDir = parameterDict['outputMorphDir'] + '/' + parameterDict['morphology'][:-4] + '/KMC'
    jobsList = [carrierList[i:i + (int(np.ceil(len(carrierList) / len(procIDs))))] for i in range(0, len(carrierList), int(np.ceil(len(carrierList)/float(len(procIDs)))))]
    print("Writing job pickles for each CPU...")
    runningJobs = []
    for procID, jobs in enumerate(jobsList):
        pickleName = outputDir + '/KMCData_%02d.pickle' % (procID)
        with open(pickleName, 'wb+') as pickleFile:
            pickle.dump(jobs, pickleFile)
        print("KMC jobs for procID", procID, "written to KMCData_%02d.pickle" % (procID))
        # Open the required processes to execute the KMC jobs
        #print('python ' + os.getcwd() + '/code/singleCoreRunMobKMC.py ' + outputDir + ' ' + str(procID) + ' &')
        runningJobs.append(sp.Popen(['python', SINGLE_RUN_MOBKMC_FILE, outputDir, str(procID)]))
    # Wait for all jobs to complete
    [p.wait() for p in runningJobs]
    # Now combine all of the pickle files into one:
    print("All KMC jobs completed!")
    if parameterDict['combineKMCResults'] is True:
        print("Combining outputs...")
        combinedData = {}
        for procID, jobs in enumerate(jobsList):
            fileName = outputDir + '/KMCResults_%02d.pickle' % (procID)
            # The pickle was repeatedly dumped to, in order to save time.
            # Each dump stream is self-contained, so iteratively unpickle to add the new data.
            with open(fileName, 'rb') as pickleFile:
                pickledData = pickle.load(pickleFile)
                for key, val in pickledData.items():
                    if key not in combinedData:
                        combinedData[key] = val
                    else:
                        combinedData[key] += val
        # Write out the combined data
        with open(outputDir + '/KMCResults.pickle', 'wb+') as pickleFile:
            pickle.dump(combinedData, pickleFile)
        print("Complete data written to", outputDir + "/KMCResults.pickle.")
        print("Cleaning up...")
        # Delete any unneeded files
        for fileName in glob.glob(outputDir + '/KMCResults_*'):
            os.remove(fileName)
        for fileName in glob.glob(outputDir + '/KMCslot_*'):
            os.remove(fileName)
    for fileName in glob.glob(outputDir + '/KMCData*'):
        os.remove(fileName)
    return AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList

if __name__ == "__main__":
    try:
        pickleFile = sys.argv[1]
    except:
        print("Please specify the pickle file to load to continue the pipeline from this point.")
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = helperFunctions.loadPickle(pickleFile)
    execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList)
