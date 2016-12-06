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


plottingSubroutines = False
elementaryCharge = 1.60217657E-19 # C
kB = 1.3806488E-23 # m^{2} kg s^{-2} K^{-1}
hbar = 1.05457173E-34 # m^{2} kg s^{-1}

# Flag for outputting all hop times to a CSV (slow!)
debug = False
###


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
                deltaPosition[axis] -= self.simDims[axis][1] - self.simDims[axis][0]
                self.image[axis] += 1
            while deltaPosition[axis] < - halfBoxLength:
                # Crossed over a negative boundary, decrement image by 1
                deltaPosition[axis] += self.simDims[axis][1] - self.simDims[axis][0]
                self.image[axis] -= 1
        # Carrier image now sorted, so update its current position
        self.currentChromophore = destinationChromophore
        # Increment the simulation time
        self.currentTime += hopTime
        # Increment the hop counter
        self.noHops += 1
        # Now update the sparse history matrix
        if self.carrierHistoryMatrix is not None:
            self.carrierHistoryMatrix[initialID, destinationID] += 1


def execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList):
    # Determine the maximum simulation times based on the parameter dictionary
    simulationTimes = parameterDict['simulationTimes']
    carrierList = []
    # For each specifed carrier
    for carrierNo in range(parameterDict['numberOfCarriersPerSimulationTime']):
        print "\rCreating carrier number", str(carrierNo + 1), "of", str(parameterDict['numberOfCarriersPerSimulationTime']) + "...",
        # For each specified simulation time (do the loops this way round to even out the job lists when we split all the KMC jobs over the CPUs)
        for lifetime in simulationTimes:
            # Find an initial position (for now, just pick a chromophore at random)
            startChromoID = R.randint(0, len(chromophoreList) - 1)
            carrierList.append(carrier(chromophoreList, parameterDict, startChromoID, lifetime, carrierNo, AAMorphologyDict))
    print ""
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
    return AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList





#    # Create pickle file containing the jobs sorted by ProcID to be picked up by singleCoreRunKMC.py
#    pickleName = outputDir + '/carrierJobs.pickle'
#    print "Writing job pickle for each CPU..."
#    with open(pickleName, 'w+') as pickleFile:
#        pickle.dump(jobsList, pickleFile)
#    print "KMC jobs list written to", pickleName
#    if len(jobsList) <= len(procIDs):
#        procIDs = procIDs[:len(jobsList)]
#    runningJobs = []
#    # Open the required processes to execute the KMC jobs
#    for CPURank in procIDs:
#        print 'python ' + os.getcwd() + '/code/singleCoreRunKMC.py ' + outputDir + ' ' + str(CPURank) + ' &'
#        runningJobs.append(sp.Popen(['python', str(os.getcwd()) + '/code/singleCoreRunKMC.py', outputDir, str(CPURank)]))
#    # Wait for all jobs to complete
#    [p.wait() for p in runningJobs]
#    # Delete the job pickle
#    os.system('rm ' + pickleName)
#    return AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList


if __name__ == "__main__":
    try:
        pickleFile = sys.argv[1]
    except:
        print "Please specify the pickle file to load to continue the pipeline from this point."
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList = helperFunctions.loadPickle(pickleFile)
    execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList)
