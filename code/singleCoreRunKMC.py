import sys
import signal
import traceback
import os
import numpy as np
import helperFunctions
import time as T
import random as R
from scipy.sparse import lil_matrix
import subprocess as sp
import pickle


elementaryCharge = 1.60217657E-19 # C
kB = 1.3806488E-23 # m^{2} kg s^{-2} K^{-1}
hbar = 1.05457173E-34 # m^{2} kg s^{-1}


class carrier:
    def __init__(self, chromophoreList, parameterDict, chromoID, lifetime, carrierNo, AAMorphologyDict):
        if parameterDict['recordCarrierHistory'] is True:
            self.carrierHistoryMatrix = lil_matrix((len(chromophoreList), len(chromophoreList)), dtype = int)
        else:
            self.carrierHistoryMatrix = None
        self.ID = carrierNo
        self.image = [0, 0, 0]
        self.initialChromophore = chromophoreList[chromoID]
        self.currentChromophore = chromophoreList[chromoID]
        if parameterDict['hopLimit'] == 0:
            self.hopLimit = None
        else:
            self.hopLimit = parameterDict['hopLimit']
        self.T = parameterDict['systemTemperature']
        self.lifetime = lifetime
        self.currentTime = 0.0
        self.lambdaij = parameterDict['reorganisationEnergy']
        self.noHops = 0
        self.simDims = [[-AAMorphologyDict['lx'] / 2.0, AAMorphologyDict['lx'] / 2.0], [-AAMorphologyDict['ly'] / 2.0, AAMorphologyDict['ly'] / 2.0], [-AAMorphologyDict['lz'] / 2.0, AAMorphologyDict['lz'] / 2.0]]
        self.displacement = None

    def calculateHopRate(self, lambdaij, Tij, deltaEij):
        # Semiclassical Marcus Hopping Rate Equation
        kij = ((2 * np.pi) / hbar) * (Tij ** 2) * np.sqrt(1.0 / (4 * lambdaij * np.pi * kB * self.T)) * np.exp(-((deltaEij + lambdaij)**2) / (4 * lambdaij * kB * self.T))
        return kij

    def determineHopTime(self, rate):
        # Use the KMC algorithm to determine the wait time to this hop
        if rate != 0:
            while True:
                x = R.random()
                # Ensure that we don't get exactly 0.0 or 1.0, which would break our logarithm
                if (x != 0.0) and (x != 1.0):
                    break
            tau = - np.log(x) / rate
        else:
            # If rate == 0, then make the hopping time extremely long
            tau = 1E99
        return tau

    def calculateHop(self, chromophoreList):
        # Terminate if the next hop would be more than the termination limit
        if self.hopLimit is not None:
            if self.noHops + 1 > self.hopLimit:
                return 1
        # Determine the hop times to all possible neighbours
        hopTimes = []
        # Obtain the reorganisation energy in J (from eV in the parameter file)
        for neighbourIndex, transferIntegral in enumerate(self.currentChromophore.neighboursTI):
            deltaEij = self.currentChromophore.neighboursDeltaE[neighbourIndex]
            # All of the energies are in eV currently, so convert them to J
            hopRate = self.calculateHopRate(self.lambdaij * elementaryCharge, transferIntegral * elementaryCharge, deltaEij * elementaryCharge)
            hopTime = self.determineHopTime(hopRate)
            # Keep track of the chromophoreID and the corresponding tau
            hopTimes.append([self.currentChromophore.neighbours[neighbourIndex][0], hopTime])
        # Sort by ascending hop time
        hopTimes.sort(key = lambda x:x[1])
        # Take the quickest hop
        destinationChromophore = chromophoreList[hopTimes[0][0]]
        # As long as we're not limiting by the number of hops:
        if self.hopLimit is None:
            # Ensure that the next hop does not put the carrier over its lifetime
            if (self.currentTime + hopTimes[0][1]) > self.lifetime:
                # Send the termination signal to singleCoreRunKMC.py
                return 1
        # Move the carrier and send the contiuation signal to singleCoreRunKMC.py
        self.performHop(destinationChromophore, hopTimes[0][1])
        return 0

    def performHop(self, destinationChromophore, hopTime):
        initialID = self.currentChromophore.ID
        destinationID = destinationChromophore.ID
        initialPosition = self.currentChromophore.posn
        destinationPosition = destinationChromophore.posn
        deltaPosition = destinationPosition - initialPosition
        for axis in range(3):
            halfBoxLength = (self.simDims[axis][1] - self.simDims[axis][0]) / 2.0
            while deltaPosition[axis] > halfBoxLength:
                # Crossed over a positive boundary, increment image by 1
                deltaPosition[axis] -= halfBoxLength * 2.0
                self.image[axis] -= 1
            while deltaPosition[axis] < - halfBoxLength:
                # Crossed over a negative boundary, decrement image by 1
                deltaPosition[axis] += halfBoxLength * 2.0
                self.image[axis] += 1
        # Carrier image now sorted, so update its current position
        self.currentChromophore = destinationChromophore
        # Increment the simulation time
        self.currentTime += hopTime
        # Increment the hop counter
        self.noHops += 1
        # Now update the sparse history matrix
        if self.carrierHistoryMatrix is not None:
            self.carrierHistoryMatrix[initialID, destinationID] += 1


class saveCarrier:
    def __init__(self, **kwargs):#ID, image, initialPosn, finalPosn, lifetime, currentTime, noHops, displacement):
        for key, value in kwargs.items():
            self.__dict__[key] = value


class terminationSignal:
    killSent = False
    def __init__(self):
        signal.signal(signal.SIGINT, self.catchKill)
        signal.signal(signal.SIGTERM, self.catchKill)

    def catchKill(self, signum, frame):
        self.killSent = True

class Terminate(Exception):
    '''This class is raised to terminate a KMC simulation'''
    def __init__(self, string):
        self.string = string

    def __str__(self):
        return self.string

def savePickle(saveData, savePickleName):
    with open(savePickleName, 'w+') as pickleFile:
        pickle.dump(saveData, pickleFile)
    helperFunctions.writeToFile(logFile, ['Pickle file saved successfully as ' + savePickleName + '!'])


def calculateDisplacement(initialPosition, finalPosition, finalImage, simDims):
    displacement = [0.0, 0.0, 0.0]
    for axis in range(3):
        displacement[axis] = (finalPosition[axis] - initialPosition[axis]) + (finalImage[axis] * (simDims[axis][1] - simDims[axis][0]))
    return np.linalg.norm(np.array(displacement))


def initialiseSaveData(nChromos, seed):
    return {'seed': seed, 'ID': [], 'image': [], 'lifetime': [], 'currentTime': [], 'noHops': [], 'displacement': [], 'carrierHistoryMatrix': lil_matrix((nChromos, nChromos), dtype = int), 'initialPosition': [], 'finalPosition': []}


if __name__ == '__main__':
    KMCDirectory = sys.argv[1]
    CPURank = int(sys.argv[2])
    overwrite = False
    try:
        overwrite = bool(sys.argv[3])
    except:
        pass
    # Load `jobsToRun' which is a list of tuples that contains the [carrier.ID, carrier.lifetime]
    pickleFileName = KMCDirectory + '/KMCData_%02d.pickle' % (CPURank)
    with open(pickleFileName, 'r') as pickleFile:
        jobsToRun = pickle.load(pickleFile)
    logFile = KMCDirectory + '/KMClog_' + str(CPURank) + '.log'
    # Reset the log file
    with open(logFile, 'w+') as logFileHandle:
        pass
    helperFunctions.writeToFile(logFile, ['Found ' + str(len(jobsToRun)) + ' jobs to run.'])
    # Set the affinities for this current process to make sure it's maximising available CPU usage
    currentPID = os.getpid()
    try:
        affinityJob = sp.Popen(['taskset', '-pc', str(CPURank), str(currentPID)], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
        # helperFunctions.writeToFile(logFile, affinityJob[0].split('\n')) #stdOut for affinity set
        # helperFunctions.writeToFile(logFile, affinityJob[1].split('\n')) #stdErr for affinity set
    except OSError:
        helperFunctions.writeToFile(logFile, ["Taskset command not found, skipping setting of processor affinity..."])
    # Now load the main morphology pickle (used as a workaround to obtain the chromophoreList without having to save it in each carrier [very memory inefficient!])
    pickleDir = KMCDirectory.replace('/KMC', '/code')
    for fileName in os.listdir(pickleDir):
        if 'pickle' in fileName:
            mainMorphologyPickleName = pickleDir + '/' + fileName
    helperFunctions.writeToFile(logFile, ['Found main morphology pickle file at ' + mainMorphologyPickleName + '! Loading data...'])
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = helperFunctions.loadPickle(mainMorphologyPickleName)
    helperFunctions.writeToFile(logFile, ['Main morphology pickle loaded!'])
    # Attempt to catch a kill signal to ensure that we save the pickle before termination
    killer = terminationSignal()
    seed = R.randint(0, sys.maxsize)
    R.seed(seed)
    # Save the pickle as a list of `saveCarrier' instances that contain the bare minimum
    saveData = initialiseSaveData(len(chromophoreList), seed)
    t0 = T.time()
    saveTime = T.time()
    saveSlot = 'slot1'
    try:
        for jobNumber, [carrierNo, lifetime] in enumerate(jobsToRun):
            t1 = T.time()
            # Find a random position to start the carrier in
            startChromoID = R.randint(0, len(chromophoreList) - 1)
            # Create the carrier instance
            thisCarrier = carrier(chromophoreList, parameterDict, startChromoID, lifetime, carrierNo, AAMorphologyDict)
            terminateSimulation = False
            while terminateSimulation is False:
                terminateSimulation = bool(thisCarrier.calculateHop(chromophoreList))
                if killer.killSent is True:
                    raise Terminate('Kill command sent, terminating KMC simulation...')
            # Now the carrier has finished hopping, let's calculate its vitals
            initialPosition = thisCarrier.initialChromophore.posn
            finalPosition = thisCarrier.currentChromophore.posn
            finalImage = thisCarrier.image
            simDims = thisCarrier.simDims
            thisCarrier.displacement = calculateDisplacement(initialPosition, finalPosition, finalImage, simDims)
            # Now the calculations are completed, create a barebones class containing the save data
            importantData = ['ID', 'image', 'lifetime', 'currentTime', 'noHops', 'displacement']
            for name in importantData:
                saveData[name].append(thisCarrier.__dict__[name])
            # Update the carrierHistoryMatrix
            if parameterDict['recordCarrierHistory'] is True:
                saveData['carrierHistoryMatrix'] += thisCarrier.carrierHistoryMatrix
            else:
                saveData['carrierHistoryMatrix'] = None
            # Then add in the initial and final positions
            saveData['initialPosition'].append(initialPosition)
            saveData['finalPosition'].append(finalPosition)
            #saveData.append(saveCarrier(**dataToSave))
            t2 = T.time()
            elapsedTime = float(t2) - float(t1)
            if elapsedTime < 60:
                timeunits = 'seconds.'
            elif elapsedTime < 3600:
                elapsedTime /= 60.0
                timeunits = 'minutes.'
            elif elapsedTime < 86400:
                elapsedTime /= 3600.0
                timeunits = 'hours.'
            else:
                elapsedTime /= 86400.0
                timeunits = 'days.'
            elapsedTime = '%.1f' % (float(elapsedTime))
            helperFunctions.writeToFile(logFile, ['Carrier hopped ' + str(thisCarrier.noHops) + ' times, over ' + str(thisCarrier.currentTime) + ' seconds, into image ' + str(thisCarrier.image) + ', for a displacement of ' + str(thisCarrier.displacement) + ', in ' + str(elapsedTime) + ' wall-clock ' + str(timeunits)])
            # Save the pickle file every hour
            if (t2 - saveTime) > 3600:
                print("Completed", jobNumber, "of", len(jobsToRun), "jobs. Making checkpoint at %3d%%" % (np.round(jobNumber + 1 / float(len(jobsToRun)) * 100)))
                helperFunctions.writeToFile(logFile, ['Completed ' + str(jobNumber) + ' jobs. Making checkpoint at %3d%%' % (np.round(jobNumber / float(len(jobsToRun)) * 100))])
                savePickle(saveData, pickleFileName.replace('Data', saveSlot + 'Results'))
                if saveSlot == 'slot1':
                    saveSlot = 'slot2'
                elif saveSlot == 'slot2':
                    saveSlot = 'slot1'
                saveTime = T.time()
    except Exception as errorMessage:
        print(traceback.format_exc())
        print("Saving the pickle file cleanly before termination...")
        helperFunctions.writeToFile(logFile, [str(errorMessage)])
        helperFunctions.writeToFile(logFile, ['Saving the pickle file cleanly before termination...'])
        savePickle(saveData, pickleFileName.replace('Data', 'TerminatedResults'))
        print("Pickle saved! Exitting Python...")
        exit()
    t3 = T.time()
    elapsedTime = float(t3) - float(t0)
    if elapsedTime < 60:
        timeunits = 'seconds.'
    elif elapsedTime < 3600:
        elapsedTime /= 60.0
        timeunits = 'minutes.'
    elif elapsedTime < 86400:
        elapsedTime /= 3600.0
        timeunits = 'hours.'
    else:
        elapsedTime /= 86400.0
        timeunits = 'days.'
    elapsedTime = '%.1f' % (float(elapsedTime))
    helperFunctions.writeToFile(logFile, ['All jobs completed in ' + elapsedTime + ' ' + timeunits])
    helperFunctions.writeToFile(logFile, ['Saving the pickle file cleanly before termination...'])
    savePickle(saveData, pickleFileName.replace('Data', 'Results'))
    helperFunctions.writeToFile(logFile, ['Exiting normally...'])
