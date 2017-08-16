import os
import glob
import sys
import numpy as np
import random as R
from scipy.sparse import lil_matrix
import pickle
import subprocess as sp
import helperFunctions


class morphologyMoiety:
    def __init__(self, molMorphName, parameterDict):
        chromophoreListLocation = parameterDict['outputMorphDir'] + '/' + molMorphName + '/code/' + molMorphName + '.pickle'
        self.AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, self.parameterDict, self.chromophoreList = helperFunctions.loadPickle(chromophoreListLocation)
        self.carrierType = self.getCarrierType()
        # Now add the occupation data to the chromophoreLists so that we can prevent double occupation in the simulations
        # The occupied property is a list that contains the device moiety coordinates where the chromophore is occupied.
        for index, chromophore in enumerate(self.chromophoreList):
            chromophore.occupied = []

    def getCarrierType(self):
        speciesPresent = []
        for chromophore in self.chromophoreList:
            speciesPresent.append(chromophore.species)
        if len(set(speciesPresent)) == 1:
            if speciesPresent[0] == 'Donor':
                return 'hole'
            elif speciesPresent[0] == 'Acceptor':
                return 'electron'
            else:
                print("Error in chromophore:")
                for key, val in chromophore.__dict__.items():
                    print(key, "=", val)
                raise SystemError("Chromophore species is neither Donor nor Acceptor")
        else:
            return 'both'


class chromophoreDataContainer:
    # A helper class that contains all of the chromophore data for ease of access from anywhere
    def __init__(self, deviceArray, moietyDictionary, wrapxy):
        self.deviceArray = deviceArray
        self.moietyDictionary = moietyDictionary
        self.wrapxy = wrapxy

    def returnChromophoreList(self, devicePosition):
        deviceMoietyType = self.deviceArray[tuple(devicePosition)]
        return self.moietyDictionary[deviceMoietyType].chromophoreList

    def returnSpecificChromophore(self, devicePosition, chromoID):
        deviceMoietyType = self.deviceArray[tuple(devicePosition)]
        return self.moietyDictionary[deviceMoietyType].chromophoreList[chromoID]

    def returnRandomChromophore(self, devicePosition):
        deviceMoietyType = self.deviceArray[tuple(devicePosition)]
        return R.choice(self.moietyDictionary[deviceMoietyType].chromophoreList)

    def returnClosestChromophoreToPosition(self, devicePosition, desiredPosition):
        closestChromoID = None
        closestDistance = 1E99

        # Check that there is an eligible device position that exists at these coordinates (i.e. there is a hop taking place within the active layer)
        # Find out which axis is out of index
        for axisNo, val in enumerate(devicePosition):
            if val >= self.deviceArray.shape[axisNo]:
                if axisNo == 2:
                    # Leaving out of the top of the device
                    return 'Top'
                elif self.wrapxy:
                    # Bring it in on the reverse side
                    devicePosition[axisNo] = 0
                    #print("Axis", axisNo, "out of bounds, setting new devicePosition to", devicePosition)
                else:
                    return 'Out of Bounds'
            if val < 0:
                if axisNo == 2:
                    # Leaving out of the bottom of the device
                    return 'Bottom'
                elif self.wrapxy:
                    # Bring it in on the reverse side
                    devicePosition[axisNo] = self.deviceArray.shape[axisNo] - 1
                    #print("Axis", axisNo, "out of bounds, setting new devicePosition to", devicePosition)
                else:
                    return 'Out of Bounds'
        deviceMoietyType = self.deviceArray[tuple(devicePosition)]
        # NOTE For efficiency, this section has been rewritten
        # The following uncommented code represents a 10x speedup over the commented code
        #for chromophore in self.moietyDictionary[deviceMoietyType].chromophoreList:
        #    separation = np.linalg.norm(np.array(desiredPosition) - np.array(chromophore.posn))
        #    if separation < closestDistance:
        #        closestDistance = separation
        #        closestChromoID = chromophore.ID
        #        continue
        positions = np.array([chromo.posn for chromo in self.moietyDictionary[deviceMoietyType].chromophoreList])
        distances = np.sqrt(np.sum((positions - np.array(desiredPosition))**2, axis=1))
        closestChromoID = np.argmin(distances)
        return self.moietyDictionary[deviceMoietyType].chromophoreList[closestChromoID]


class morphologyDataContainer:
    # A helper class that contains all of the chromophore data for ease of access from anywhere
    def __init__(self, deviceArray, moietyDictionary):
        self.deviceArray = deviceArray
        self.moietyDictionary = moietyDictionary

    def returnAAMorphology(self, devicePosition):
        deviceMoietyType = self.deviceArray[tuple(devicePosition)]
        return self.moietyDictionary[deviceMoietyType].AAMorphologyDict

    def returnDeviceMoietyType(self, devicePosition):
        return self.deviceArray[tuple(devicePosition)]


def loadDeviceMorphology(parameterDict):
    deviceDir = parameterDict['inputDeviceDir'] + '/' + parameterDict['deviceMorphology']
    ySlices = os.listdir(deviceDir)
    # Initialize the array of the correct size (assumes cubic morphology)
    deviceArray = np.zeros([len(ySlices)]*3, dtype=int)
    for yVal, fileName in enumerate(ySlices):
        # Load the ySlice as-presented in the input files
        ySlice = np.loadtxt(deviceDir + '/' + fileName, dtype=int)
        # The z-origin is at the top, and we need it at the bottom, so turn the array upside down
        ySlice = np.flipud(ySlice)
        # Now populate the array
        for zVal, zRow in enumerate(ySlice):
            for xVal, datum in enumerate(zRow):
                deviceArray[xVal, yVal, zVal] = datum
    moietyDictionary = {}
    for moietyID in np.unique(deviceArray):
        moietyDictionary[moietyID] = morphologyMoiety(parameterDict['deviceComponents'][moietyID], parameterDict)
    return deviceArray, moietyDictionary


def execute(parameterDict):
    # Set the random seed now for all the child processes
    R.seed(3232)

    # First job will be to load in the device morphology, when I work out what format I want it to be.
    deviceArray, moietyDictionary = loadDeviceMorphology(parameterDict)

    # Initialise the helperClass to obtain all of the chromophoreData required, allowing it be accessed globally
    chromophoreData = chromophoreDataContainer(deviceArray, moietyDictionary, parameterDict['wrapDeviceXY'])
    morphologyData = morphologyDataContainer(deviceArray, moietyDictionary)

    # Write these classes out to a pickle file so that they can be loaded by the child processes later
    toPickle = [deviceArray, chromophoreData, morphologyData, parameterDict]
    saveDirectory = parameterDict['outputDeviceDir'] + '/' + parameterDict['deviceMorphology'] + '/code'
    if parameterDict['overwriteCurrentData'] is True:
        with open(saveDirectory + '/deviceData.pickle', 'wb+') as pickleFile:
            pickle.dump(toPickle, pickleFile)
    voltages = []
    for V in parameterDict['voltageSweep']:
        voltages.append(V)
    procIDs = parameterDict['procIDs']
    jobsList = [voltages[i:i + (int(np.ceil(len(voltages) / len(procIDs))))] for i in range(0, len(voltages), int(np.ceil(len(voltages)/float(len(procIDs)))))]
    runningJobs = []
    outputDir = parameterDict['outputDeviceDir'] + '/' + parameterDict['deviceMorphology'] + '/KMC'
    print("Writing job pickles for each CPU...")
    for procID, jobs in enumerate(jobsList):
        pickleName = outputDir + '/KMCData_%02d.pickle' % (procID)
        with open(pickleName, 'wb+') as pickleFile:
            pickle.dump(jobs, pickleFile)
        print("KMC jobs for procID", procID, "written to KMCData_%02d.pickle" % (procID))
        # Open the required processes to execute the KMC jobs
        # Random seeding is a little weird here. If we don't generate a random seed in the child process, it will just use the system time. So, we generate a seed here to get the same random number stream each time, and then feed the child process a new seed from the random number stream. This way, we ensure that each child process has a different random number stream to the other processes, but it's the same stream every time we run the program.
        childSeed = R.randint(0, 2**32)
        # NOTE DEBUG USING cProfile
        #print('RUNNING CPROFILE...')
        #print('python -m cProfile -o deviceSim.cprof ' + os.getcwd() + '/code/singleCoreRunDeviceKMC.py', outputDir + ' ' + str(procID) + ' ' + str(childSeed) + ' &')
        #runningJobs.append(sp.Popen(['python -m cProfile -o deviceSim.cprof ' + str(os.getcwd()) + '/code/singleCoreRunDeviceKMC.py', outputDir, str(procID), str(childSeed)]))
        # Previous run command:
        print('python ' + os.getcwd() + '/code/singleCoreRunDeviceKMC.py' + outputDir + ' ' + str(procID) + ' ' + str(childSeed) + ' &')
        runningJobs.append(sp.Popen(['python', str(os.getcwd()) + '/code/singleCoreRunDeviceKMC.py', outputDir, str(procID), str(childSeed)]))
    # Wait for all jobs to complete
    [p.wait() for p in runningJobs]
    print("All KMC jobs completed!")
    # Combine results if required.

if __name__ == "__main__":
    try:
        pickleFile = sys.argv[1]
    except:
        print("Please specify the pickle file to load to continue the pipeline from this point.")
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = helperFunctions.loadPickle(pickleFile)
    execute(parameterDict)
