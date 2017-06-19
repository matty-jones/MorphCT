import os
import sys
import copy
import random as R
import numpy as np
import heapq
import matplotlib.pyplot as plt
try:
    import mpl_toolkits.mplot3d as p3
except ImportError:
    print("Could not import 3D plotting engine, calling the plotDeviceComponents function will result in an error!")

import helperFunctions


# Physical Constants
elementaryCharge = 1.60217657E-19  # C
kB = 1.3806488E-23                 # m^{2} kg s^{-2} K^{-1}
hbar = 1.05457173E-34              # m^{2} kg s^{-1}
lightSpeed = 299792458             # ms^{-1}

# Global Variables
globalTime = 0
chromophoreData = []
morphologyData = []
numberOfExtractions = 0


class morphologyMoiety:
    def __init__(self, moietyTypeNumber, molMorphName, parameterDict):
        self.typeNo = moietyTypeNumber
        chromophoreListLocation = parameterDict['outputMorphDir'] + '/' + molMorphName + '/code/' + molMorphName + '.pickle'
        self.AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, self.parameterDict, self.chromophoreList = helperFunctions.loadPickle(chromophoreListLocation)
        self.carrierType = self.getCarrierType()

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
    def __init__(self, deviceArray, moietyDictionary):
        self.deviceArray = deviceArray
        self.moietyDictionary = moietyDictionary

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
                    print("Leaving out of top of device")
                    return 'Top'
                else:
                    # Bring it in on the reverse side
                    devicePosition[axisNo] = 0
                    #print("Axis", axisNo, "out of bounds, setting new devicePosition to", devicePosition)
            if val < 0:
                if axisNo == 2:
                    # Leaving out of the bottom of the device
                    print("Leaving out of bottom of device")
                    return 'Bottom'
                else:
                    # Bring it in on the reverse side
                    devicePosition[axisNo] = self.deviceArray.shape[axisNo] - 1
                    #print("Axis", axisNo, "out of bounds, setting new devicePosition to", devicePosition)
        deviceMoietyType = self.deviceArray[tuple(devicePosition)]
        for chromophore in self.moietyDictionary[deviceMoietyType].chromophoreList:
            separation = helperFunctions.calculateSeparation(desiredPosition, chromophore.posn)
            if separation < closestDistance:
                closestDistance = separation
                closestChromoID = chromophore.ID
                continue
        return self.moietyDictionary[deviceMoietyType].chromophoreList[closestChromoID]


class morphologyDataContainer:
    # A helper class that contains all of the chromophore data for ease of access from anywhere
    def __init__(self, deviceArray, moietyDictionary):
        self.deviceArray = deviceArray
        self.moietyDictionary = moietyDictionary

    def returnAAMorphology(self, devicePosition):
        deviceMoietyType = self.deviceArray[tuple(devicePosition)]
        return self.moietyDictionary[deviceMoietyType].AAMorphologyDict


class exciton:
    def __init__(self, index, globalTime, initialPosnDevice, parameterDict):
        # Initialise the important variables to be used later
        self.ID = index
        self.creationTime = globalTime
        self.removedTime = None
        self.T = parameterDict['systemTemperature']
        self.lifetimeParameter = parameterDict['excitonLifetime']
        self.recombinationTime = -np.log(R.random()) * self.lifetimeParameter
        self.rF = parameterDict['forsterRadius']
        self.numberOfHops = 0
        # NOTE For now, we'll just inject it randomly somewhere in the system
        self.initialDevicePosn = initialPosnDevice
        self.currentDevicePosn = copy.deepcopy(initialPosnDevice)
        self.initialChromophore = chromophoreData.returnRandomChromophore(initialPosnDevice)
        self.currentChromophore = copy.deepcopy(self.initialChromophore)
        self.calculateBehaviour()

    def calculateBehaviour(self):
        if globalTime >= self.creationTime + self.recombinationTime:
            # Exciton has run out of time and is now recombining, so set its next hopTime to None
            [self.destinationChromophore, self.hopTime, self.destinationImage] = [None, None, None]
            self.removedTime = globalTime
            return
        # Set a flag to indicate whether exciton can dissociate or not
        self.canDissociate = self.checkDissociation()
        # Dissociate the exciton immediately after creation if it would be created at an interface
        if self.canDissociate is True:
            # Determine all potential dissociation options, randomly select one, plop a hole on the donor and electron on the acceptor, then remove this exciton from the system by updating its removedTime
            if self.currentChromophore.species == 'Donor':
                self.holeChromophore = self.currentChromophore
                # The electron chromophore is a randomly selected chromophore of the opposing type that is in range of the current one
                electronChromophoreID = R.choice(self.currentChromophore.dissociationNeighbours)[0]
                self.electronChromophore = chromophoreData.returnSpecificChromophore(self.currentDevicePosn, electronChromophoreID)
            elif self.currentChromophore.species == 'Acceptor':
                self.electronChromophore = self.currentChromophore
                # The hole chromophore is a randomly selected chromophore of the opposing type that is in range of the current one
                holeChromophoreID = R.choice(self.currentChromophore.dissociationNeighbours)[0]
                self.holeChromophore = chromophoreData.returnSpecificChromophore(self.currentDevicePosn, holeChromophoreID)
            # Notify execute() that this exciton should not be queued up again by setting self.hopTime == None
            [self.destinationChromophore, self.hopTime, self.destinationImage] = [None, None, None]
            self.removeTime = globalTime
        # Calculate the fastest hop from the current chromophore
        try:
            [self.destinationChromophore, self.hopTime, self.destinationImage] = self.calculateHop()
        except ValueError:
            # To remove the exciton from the system (or, more acurately, to prevent the code from queueing up the exciton again because it has already been popped from the main KMC queue), self.hopTime needs to be None
            [self.destinationChromophore, self.hopTime, self.destinationImage] = [None, None, None]
            self.removedTime = globalTime

    def checkDissociation(self):
        # Return True if there are neighbours of the opposing electronic species present, otherwise return False
        if len(self.currentChromophore.dissociationNeighbours) > 0:
            return True
        return False

    def calculateHopRate(self, rij, deltaEij):
        # Foerster Transport Hopping Rate Equation
        if deltaEij <= 0:
            boltzmannFactor = 1
        else:
            boltzmannFactor = np.exp(-(elementaryCharge * deltaEij)/(kB * self.T))

        kFRET = (1/self.lifetimeParameter) * (self.rF / rij)**6 * boltzmannFactor
        return kFRET

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

    def calculateHop(self):
        # First thing's first, if the current global time is > creationTime + lifeTime then this exciton recombined before the hop so we need to remove it from the system
        if globalTime > (self.creationTime + self.recombinationTime):
            return []

        # Get the chromophoreList so that we can work out which chromophores to hop to
        chromophoreList = chromophoreData.returnChromophoreList(self.currentDevicePosn)
        # Determine the hop times to all possible neighbours
        hopTimes = []
        for neighbourIndex, transferIntegral in enumerate(self.currentChromophore.neighboursTI):
            # Ignore any hops with a NoneType transfer integral (usually due to an ORCA error)
            if transferIntegral is None:
                continue
            # For the hop, we need to know the change in Eij for the boltzmann factor, as well as the distance rij being hopped
            # Durham code has a prefactor of 1.414 here, not sure why
            deltaEij = self.currentChromophore.neighboursDeltaE[neighbourIndex]
            # The current posn inside the wrapped morphology is easy
            currentChromoPosn = self.currentChromophore.posn
            # The hop destination in wrapped morphology is slightly more complicated
            neighbourPosn = chromophoreList[self.currentChromophore.neighbours[neighbourIndex][0]].posn
            neighbourID = chromophoreList[self.currentChromophore.neighbours[neighbourIndex][0]].ID
            neighbourUnwrappedPosn = chromophoreList[self.currentChromophore.neighbours[neighbourIndex][0]].unwrappedPosn
            neighbourImage = chromophoreList[self.currentChromophore.neighbours[neighbourIndex][0]].image
            # Need to determine the relative image between the two (recorded in obtainChromophores). We will also need this to check if the exciton is hopping out of this morphology cell and into an adjacent one
            neighbourRelativeImage = self.currentChromophore.neighbours[neighbourIndex][1]
            #if neighbourRelativeImage != [0, 0, 0] and neighbourID == 1063:
            #    print("Current =", self.currentChromophore.ID, currentChromoPosn)
            #    print("Destination =", neighbourID, neighbourPosn)
            #    print(neighbourPosn, neighbourUnwrappedPosn, neighbourRelativeImage)
            #    morphologyShape = [morphologyData.returnAAMorphology(self.currentDevicePosn)[key] for key in ['lx', 'ly', 'lz']]
            #    print(neighbourPosn + [neighbourRelativeImage[axis] * morphologyShape[axis] for axis in range(3)])
            #    input("PAUSE...")
            # Morphology shape so we can work out the actual relative positions
            morphologyShape = [morphologyData.returnAAMorphology(self.currentDevicePosn)[key] for key in ['lx', 'ly', 'lz']]
            remappedPosn = neighbourPosn + [neighbourRelativeImage[axis] * morphologyShape[axis] for axis in range(3)]
            rij = helperFunctions.calculateSeparation(currentChromoPosn, remappedPosn) * 1E-10
            #if neighbourRelativeImage != [0, 0, 0]:
            #    print(neighbourID, neighbourPosn, neighbourUnwrappedPosn, neighbourImage, neighbourRelativeImage, rij)
            # Note, separations are recorded in angstroems, so spin this down to metres.
            # Additionally, all of the energies are in eV currently, so convert them to J
            hopRate = self.calculateHopRate(rij, deltaEij * elementaryCharge)
            hopTime = self.determineHopTime(hopRate)
            # Keep track of the destination chromophore ID, the corresponding tau, and the relative image (to see if we're hopping over a boundary)
            hopTimes.append([chromophoreData.returnChromophoreList(self.currentDevicePosn)[self.currentChromophore.neighbours[neighbourIndex][0]], hopTime, neighbourRelativeImage])
        # Sort by ascending hop time
        hopTimes.sort(key = lambda x:x[1])
        # Only want to take the quickest hop if it is NOT going to hop outside the device (top or bottom - i.e. z axis is > shape[2] or < 0)
        # Other-axis periodic hops (i.e. x and y) are permitted, and will allow the exciton to loop round again.
        while len(hopTimes) > 0:
            # Get the hop destination cell
            newLocation = list(np.array(self.currentDevicePosn) + np.array(hopTimes[0][2]))
            # If the z axis component is > shape[2] or < 0, then forbid this hop by removing it from the hopTimes list.
            if (newLocation[2] >= chromophoreData.deviceArray.shape[2]) or (newLocation[2] < 0):
                hopTimes.pop(0)
            else:
                break
        # Take the quickest hop
        if len(hopTimes) > 0:
            return hopTimes[0]
        else:
            # The exciton has not yet recombined, but there are no legal hops that can be performed so the only fate for this exciton is recombination through photoluminescence 
            return []

    def performHop(self):
        initialID = self.currentChromophore.ID
        destinationID = self.destinationChromophore.ID
        initialPosition = self.currentChromophore.posn
        destinationPosition = self.destinationChromophore.posn
        deltaPosition = destinationPosition - initialPosition
        if self.destinationImage == [0, 0, 0]:
            # Exciton is not hopping over a boundary, so can simply update its current position
            self.currentChromophore = self.destinationChromophore
        else:
            #print("Hopping into", np.array(self.currentDevicePosn) + np.array(self.destinationImage))
            # We're hopping over a boundary. Permit the hop with the already-calculated hopTime, but then find the closest chromophore to the destination position in the adjacent device cell.
            self.currentDevicePosn = list(np.array(self.currentDevicePosn) + np.array(self.destinationImage))
            newChromophore = chromophoreData.returnClosestChromophoreToPosition(self.currentDevicePosn, destinationPosition)
            if (newChromophore == 'Top') or (newChromophore == 'Bottom'):
                # This exciton is hopping out of the active layer and into the contacts.
                # Ensure it doesn't get queued up again
                [self.destinationChromophore, self.hopTime, self.destinationImage] = [None, None, None]
                self.removedTime = globalTime
                # TODO This carrier has left the simulation volume so we now need to ensure we remove it from the carrier dictionary so it's not included in any of the Coulombic calculations
            else:
                self.currentChromophore = newChromophore
        self.numberOfHops += 1
        self.canDissociate = self.checkDissociation()


class carrier:
    def __init__(self, index, globalTime, initialPosnDevice, initialChromophore, injectedFrom, parameterDict):
        self.ID = index
        self.creationTime = globalTime
        self.removedTime = None
        self.initialDevicePosn = initialPosnDevice
        self.currentDevicePosn = copy.deepcopy(initialPosnDevice)
        self.initialChromophore = initialChromophore
        self.currentChromophore = copy.deepcopy(initialChromophore)
        self.injectedFrom = injectedFrom
        self.T = parameterDict['systemTemperature']
        if self.currentChromophore.species == 'Donor':
            self.carrierType = 'Hole'
            self.lambdaij = parameterDict['reorganisationEnergyDonor']
        elif self.currentChromophore.species == 'Acceptor':
            self.carrierType = 'Electron'
            self.lambdaij = parameterDict['reorganisationEnergyAcceptor']
        self.calculateBehaviour()

    def calculateBehaviour(self):
        [self.destinationChromophore, self.hopTime, self.destinationImage] = self.calculateHop()
        try:
            [self.destinationChromophore, self.hopTime, self.destinationImage] = self.calculateHop()
        except:
            [self.destinationChromophore, self.hopTime, self.destinationImage] = [None, None, None]

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

    def calculateHop(self):
        # Determine the hop times to all possible neighbours
        hopTimes = []
        # Obtain the reorganisation energy in J (from eV in the parameter file)
        for neighbourIndex, transferIntegral in enumerate(self.currentChromophore.neighboursTI):
            # Ignore any hops with a NoneType transfer integral (usually due to an ORCA error)
            if transferIntegral is None:
                continue
            # TODO Need to include the additional deltaEij terms here
            deltaEij = self.currentChromophore.neighboursDeltaE[neighbourIndex]
            # All of the energies are in eV currently, so convert them to J
            hopRate = self.calculateHopRate(self.lambdaij * elementaryCharge, transferIntegral * elementaryCharge, deltaEij * elementaryCharge)
            hopTime = self.determineHopTime(hopRate)
            # Keep track of the chromophoreID and the corresponding tau
            hopTimes.append([self.currentChromophore.neighbours[neighbourIndex][0], hopTime])
        # Sort by ascending hop time
        hopTimes.sort(key = lambda x:x[1])
        # Take the quickest hop
        if len(hopTimes) > 0:
            return hopTimes[0]
        else:
            # We are trapped here, so just return in order to raise the calculateBehaviour exception
            return []

    def performHop(self):
        global numberOfExtractions
        initialID = self.currentChromophore.ID
        destinationID = self.destinationChromophore.ID
        initialPosition = self.currentChromophore.posn
        destinationPosition = self.destinationChromophore.posn
        deltaPosition = destinationPosition - initialPosition
        if self.destinationImage == [0, 0, 0]:
            # Carrier is not hopping over a boundary, so can simply update its current position
            self.currentChromophore = self.destinationChromophore
        else:
            # We're hopping over a boundary. Permit the hop with the already-calculated hopTime, but then find the closest chromophore to the destination position in the adjacent device cell.
            self.currentDevicePosn = list(np.array(self.currentDevicePosn) + np.array(self.destinationImage))
            newChromophore = chromophoreData.returnClosestChromophoreToPosition(self.currentDevicePosn, destinationPosn)
            if (newChromophore == 'Top') or (newChromophore == 'Bottom'):
                # This carrier is hopping out of the active layer and into the contacts.
                # Firstly, ensure that it doesn't get queued up again
                [self.destinationChromophore, self.hopTime, self.destinationImage] = [None, None, None]
                self.removedTime = globalTime
                # Secondly, work out whether this is a `correct' hop (i.e. hole hopping to anode or electron hopping to cathode) that causes photovoltaic current.
                if newChromophore == 'Top':
                    # Leaving through top (anode)
                    if self.injectedFrom != 'Anode':
                        if self.carrierType == 'Hole':
                            numberOfExtractions += 1
                        else:
                            numberOfExtractions -= 1
                    # Else (injected from anode), number of extractions doesn't change.
                else:
                    # Leaving through bottom (cathode)
                    if self.injectedFrom != 'Cathode':
                        if self.carrierType == 'Electron':
                            numberOfExtractions += 1
                        else:
                            numberOfExtractions -= 1
                    # Else (injected from cathode), number of extractions doesn't change.
            else:
                self.currentChromophore = newChromophore


def plotHopDistance(distribution):
    plt.figure()
    plt.hist(distribution)
    plt.show()
    exit()


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
        moietyDictionary[moietyID] = morphologyMoiety(moietyID, parameterDict['deviceComponents'][moietyID], parameterDict)
    return deviceArray, moietyDictionary



def plotDeviceComponents(deviceArray):
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    colour = ['r', 'g', 'b', 'y', 'k']
    for xVal in range(9):
        for yVal in range(9):
            for zVal in range(9):
                ax.scatter(xVal, yVal, zVal, zdir='z', c= colour[deviceArray[xVal, yVal, zVal]])
    plt.show()
    exit()


def calculatePhotoinjectionRate(parameterDict, deviceShape):
    # Photoinjection rate is given by the following equation. Calculations will be performed in SI always.
    rate = ((parameterDict['incidentFlux'] * 10) *                                                 # mW/cm^{2} to W/m^{2}
           (parameterDict['incidentWavelength'] / (hbar * 2 * np.pi * lightSpeed)) *               # already in m
            deviceShape[0] * deviceShape[1] * parameterDict['morphologyCellSize']**2 *             # area of device in m^{2}
           (1 - np.exp(-100 * parameterDict['absorptionCoefficient'] * deviceShape[2] * parameterDict['morphologyCellSize'])))
    return rate


def calculateDarkInjectionRate(parameterDict, deviceShape):
    # No Dark current just yet - it'll be hard to determine
    return 1E99


def decrementTime(eventQueue, eventTime):
    global globalTime
    # A function that increments the global time by whatever the time this event is, and
    # then decrements all of the remaining times in the queue.
    globalTime += eventTime
    eventQueue = [(event[0] - eventTime, event[1], event[2]) for event in eventQueue]
    return eventQueue




def plotConnections(excitonPath, chromophoreList, AADict):
    # A complicated function that shows connections between carriers in 3D that carriers prefer to hop between.
    # Connections that are frequently used are highlighted in black, whereas rarely used connections are more white.
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    for index, chromophoreID in enumerate(excitonPath[:-1]):
        coords1 = chromophoreList[chromophoreID].posn
        coords2 = chromophoreList[excitonPath[index + 1]].posn
        ax.plot([coords1[0], coords2[0]], [coords1[1], coords2[1]], [coords1[2], coords2[2]], color = 'r', linewidth = 0.5)
    simDims = [AADict['lx']/2.0, AADict['ly']/2.0, AADict['lz']/2.0]
    corners = []
    corners.append([-simDims[1], -simDims[1], -simDims[2]])
    corners.append([-simDims[0], simDims[1], -simDims[2]])
    corners.append([-simDims[0], -simDims[1], simDims[2]])
    corners.append([-simDims[0], simDims[1], simDims[2]])
    corners.append([simDims[0], -simDims[1], -simDims[2]])
    corners.append([simDims[0], simDims[1], -simDims[2]])
    corners.append([simDims[0], -simDims[1], simDims[2]])
    corners.append([simDims[0], simDims[1], simDims[2]])
    for corner1 in corners:
        for corner2 in corners:
            ax.plot([corner1[0], corner2[0]], [corner1[1], corner2[1]], [corner1[2], corner2[2]], color = 'k')
    fileName = '3d_exciton.pdf'
    plt.show()
    plt.savefig('./' + fileName)
    print("Figure saved as ./" + fileName)


def execute(parameterDict):
    # ---=== PROGRAMMER'S NOTE ===---
    # This `High Resolution' version of the code will permit chromophore-based hopping through the device, rather than just approximating the distribution. We have the proper resolution there, let's just use it!
    # ---=========================---
    global chromophoreData
    global morphologyData
    global globalTime

    R.seed(3232)

    # First job will be to load in the device morphology, when I work out what format I want it to be.
    deviceArray, moietyDictionary = loadDeviceMorphology(parameterDict)

    # Initialise the helperClass to obtain all of the chromophoreData required, allowing it be accessed globally
    chromophoreData = chromophoreDataContainer(deviceArray, moietyDictionary)
    morphologyData = morphologyDataContainer(deviceArray, moietyDictionary)

    # DEBUG
    #morphCellSizes = []
    #for moietyType in morphologyData.moietyDictionary.keys():
    #    morphCellSizes.append(morphologyData.moietyDictionary[moietyType].AAMorphologyDict['lx'])
    #print(np.average(morphCellSizes))
    #exit()

    # Need to initalise a bunch of variables
    eventQueue = []
    # Need a carrier dictionary so that we can work out energetic landscape and remove the right carriers
    carriers = {}
    photoinjectionRate = calculatePhotoinjectionRate(parameterDict, deviceArray.shape)
    darkInjectionRate = calculateDarkInjectionRate(parameterDict, deviceArray.shape)
    outputCurrentDensity = []
    numberOfPhotoinjections = 0
    excitonIndex = 0
    carrierIndex = 0
    outputCurrentConverged = False

    # Calculate the convergence characteristics
    checkConvEvery = int(parameterDict['minimumNumberOfPhotoinjections'] / 100)*5
    previousCheck = numberOfPhotoinjections
    convGlobalTime = []
    convExtractions = []
    convGradients = []



    # DEBUG
    dissExcitonDisp = []
    dissExcitonTime = []
    recExcitonDisp = []
    recExcitonTime = []
    numberOfHops = []
    numberOfDissociations = 0
    numberOfRecombinations = 0


    print("\n\n ---=== MAIN KMC LOOP START ===---\n\n")
    # Main KMC loop
    while True:
        #print("\rGlobal Time =", globalTime, "Queue =", eventQueue)
        # TODO Check whether the output current has converged, after we have had sufficient photoinjections
        # To check convergence, create a datapoint every time 5% of the parameter-file-determined photoinjections
        # have been completed. This datapoint is the numberOfExtractions as a function of time. When the required number of
        # photoinjections has been completed, start performing linear regression on the dataset and set outputCurrentConverged
        # when the numberOfExtractions saturates.
        if numberOfPhotoinjections >= (previousCheck + checkConvEvery):
            # Regardless of how many photoinjections we have, add to the convergence datasets every 5% of the minimum required
            convGlobalTime.append(globalTime)
            convExtractions.append(numberOfExtractions)
            if numberOfPhotoinjections > parameterDict['minimumNumberOfPhotoinjections']:
                # If we've gone over the minimum, then perform the convergence check to see if the output current has converged
                # TODO Check convergence, but for now, let's just see what the curve looks like
                break

        if outputCurrentConverged is True:
            break

        if len(eventQueue) == 0:
            # Either the simulation has just started, or everything just recombined/left the device
            # Therefore, queue up one photoinjection and one dark current injection to get us rolling.
            # These are the only available kinetic starting points for the simulation. Everything else stems from these.
            # Photoinjection
            photoinjectionTime = -np.log(R.random()) / photoinjectionRate
            heapq.heappush(eventQueue, (photoinjectionTime, 'photo', None))
            numberOfPhotoinjections += 1
            # Dark Injection:
            # TODO


        #print("Event Queue Before Selection =", eventQueue)

        # Now find out what the next event is
        nextEvent = heapq.heappop(eventQueue)
        # Increment the global time and decrement all of the other times in the queue
        #print("Event Queue after event has been popped", eventQueue)
        eventQueue = decrementTime(eventQueue, nextEvent[0])
        #print("Event Queue after time has been decremented", eventQueue)

        # Execute the next behaviour (being sure to recalculate rates for any new particles that have appeared)
        if nextEvent[1] == 'photo':
            # Complete the event by injecting an exciton
            # First find an injection location. For now, this will just be somewhere random in the system
            randomDevicePosition = [R.randint(0, x-1) for x in deviceArray.shape]
            print("EVENT: Photoinjection into", randomDevicePosition, "which has type", deviceArray[tuple(randomDevicePosition)])
            injectedExciton = exciton(excitonIndex, globalTime, randomDevicePosition, parameterDict)
            if (injectedExciton.canDissociate is True) or (injectedExciton.hopTime is None):
                # Injected onto either a dissociation site or a trap site
                if hoppingExciton.canDissociate is True:
                    print("EVENT: Exciton Dissociating!")
                    numberOfDissociations += 1
                    numberOfHops.append(injectedExciton.numberOfHops)
                    #injectedElectron = carrier(carrierIndex, globalTime, injectedExciton.currentDevicePosn, injectedExciton.electronChromophore, 'Exciton', parameterDict)
                    #carrier[carrierIndex] = injectedElectron
                    #carrierIndex += 1
                    #injectedHole = carrier(carrierIndex, globalTime, injectedExciton.currentDevicePosn, injectedExciton.holeChromophore, 'Exciton', parameterDict)
                    ## Add the hole to the carrier list for when we need to calc deltaE in the device
                    #carrier[carrierIndex] = injectedHole
                    #carrierIndex += 1
                    ## Now add both carriers' next hops to the KMC queue
                    #heapq.heappush(eventQueue, (injectedElectron.hopTime, 'carrierHop', injectedElectron))
                    #heapq.heappush(eventQueue, (injectedHole.hopTime, 'carrierHop', injectedHole))
                else:
                    # Injected onto a site with no connections, so this exciton will eventually die
                    print("EVENT: Exciton Recombining")
                    numberOfRecombinations += 1
                    numberOfHops.append(injectedExciton.numberOfHops)
            else:
                # Hopping permitted, push the exciton to the queue.
                # Push the exciton to the queue
                heapq.heappush(eventQueue, (injectedExciton.hopTime, 'excitonHop', injectedExciton))
                # Increment the exciton counter
                excitonIndex += 1
            # A photoinjection has just occured, so now queue up a new one
            photoinjectionTime = -np.log(R.random()) / photoinjectionRate
            print("Current Queue =", eventQueue, "Adding in:", (photoinjectionTime, 'photo', None))
            heapq.heappush(eventQueue, (photoinjectionTime, 'photo', None))
            print("Added photoinjection number", numberOfPhotoinjections)
            numberOfPhotoinjections += 1

        elif nextEvent[1] == 'dark':
            raise SystemError("No dark injection yet")

        elif nextEvent[1] == 'excitonHop':
            #print("EVENT: Exciton Hopping in device coordinates", nextEvent[2].currentDevicePosn, "(type =", deviceArray[tuple(nextEvent[2].currentDevicePosn)], ") from chromophore", nextEvent[2].currentChromophore.ID, "to chromophore", nextEvent[2].destinationChromophore.ID)
            #excitonPath.append(nextEvent[2].currentChromophore.ID)
            #excitonTime.append(globalTime - nextEvent[2].creationTime)
            #chromoList = chromophoreData.returnChromophoreList(nextEvent[2].currentDevicePosn)
            #currentChromoPosn = chromoList[nextEvent[2].currentChromophore.ID].posn
            #destinationChromoPosn = chromoList[nextEvent[2].destinationChromophore.ID].posn
            #excitonDisp.append(helperFunctions.calculateSeparation(currentChromoPosn, destinationChromoPosn))
            #if len(excitonTime) == 1000:
            #    break

            hoppingExciton = nextEvent[2]
            hoppingExciton.performHop()
            if hoppingExciton.removedTime is None:
                # As long as this exciton hasn't just been removed, recalculate its behaviour
                hoppingExciton.calculateBehaviour()
            # At this point, dissociate the exciton or remove it from the system if needed.
            if (hoppingExciton.canDissociate is True) or (hoppingExciton.hopTime is None):
                # Exciton needs to be removed. As we've already popped it from the queue, we just need to not queue it up again.
                if hoppingExciton.canDissociate is True:
                    print("EVENT: Exciton Dissociating!")
                    numberOfDissociations += 1
                    numberOfHops.append(injectedExciton.numberOfHops)
                    #input("Pause...")
                    #injectedElectron = carrier(carrierIndex, globalTime, hoppingExciton.currentDevicePosn, hoppingExciton.electronChromophore, 'Exciton', parameterDict)
                    ## Add the electron to the carrier list for when we need to calc deltaE in the device
                    #carriers[carrierIndex] = injectedElectron
                    #carrierIndex += 1
                    #injectedHole = carrier(carrierIndex, globalTime, hoppingExciton.currentDevicePosn, hoppingExciton.holeChromophore, 'Exciton', parameterDict)
                    ## Add the hole to the carrier list for when we need to calc deltaE in the device
                    #carriers[carrierIndex] = injectedHole
                    #carrierIndex += 1
                    ## Now add both carriers' next hops to the KMC queue
                    #print(eventQueue)
                    #print(injectedElectron.hopTime, injectedHole.hopTime)
                    #heapq.heappush(eventQueue, (injectedElectron.hopTime, 'carrierHop', injectedElectron))
                    #heapq.heappush(eventQueue, (injectedHole.hopTime, 'carrierHop', injectedHole))
                else:
                    print("EVENT: Exciton Recombining")
                    numberOfRecombinations += 1
                    numberOfHops.append(injectedExciton.numberOfHops)
                    #pass
                    #input("Pause...")
                # DEBUG
                # Calculate the initial position and final positions and append the excitonDisp with the separation
                initialPos = np.array(hoppingExciton.initialDevicePosn)*parameterDict['morphologyCellSize'] + (np.array(hoppingExciton.initialChromophore.posn) * 1E-10)
                finalPos = np.array(hoppingExciton.currentDevicePosn)*parameterDict['morphologyCellSize'] + (np.array(hoppingExciton.currentChromophore.posn) * 1E-10)
                if hoppingExciton.canDissociate is True:
                    dissExcitonDisp.append(helperFunctions.calculateSeparation(initialPos, finalPos) / 1E-9)
                    dissExcitonTime.append(hoppingExciton.recombinationTime)
                else:
                    recExcitonDisp.append(helperFunctions.calculateSeparation(initialPos, finalPos) / 1E-9)
                    recExcitonTime.append(hoppingExciton.recombinationTime)
                # END DEBUG
            else:
                heapq.heappush(eventQueue, (hoppingExciton.hopTime, 'excitonHop', injectedExciton))

        elif nextEvent[1] == 'carrierHop':
            raise SystemError("FOR NOW GET EXCITONS TO WORK")


            print("EVENT: Carrier Hopping")
            print(eventQueue)
            input("Pause...")
            hoppingCarrier = nextEvent[2]
            hoppingCarrier.performHop()
            hoppingCarrier.calculateBehaviour()
            if hoppingCarrier.hopTime is not None:
                heapq.heappush(eventQueue, (hoppingCarrier.hopTime, 'carrierHop', hoppingCarrier))
            # Else: The carrier has become stuck (this should only be possible if we've crossed a device cell boundary and dropped the carrier on a chromophore with no corrections).
            # The best way to deal with this is to not queue it up, but leave it in the carrier dictionary so that it still affects the energetic landscape and can be recombined with.
            # Otherwise it will just sit here for the rest of the simulation and not do anything. Not sure if I like that.
            print("Queue after carrier has hopped")
            print(eventQueue)
            input("Pause...")

        elif nextEvent[1] == 'recombine':
            raise SystemError("No recombine yet")

        else:
            print(eventQueue)
            raise SystemError("New Event is next in queue")

    print("Plotting Exciton Characteristics")
    plt.figure()
    plt.hist(dissExcitonDisp, bins=10)
    plt.title('dissExcitonDisp')
    plt.savefig('./dissExcitonDisp.pdf')
    plt.clf()

    plt.hist(dissExcitonTime, bins=10)
    plt.title('dissExcitonTime')
    plt.savefig('./dissExcitonTime.pdf')
    plt.clf()

    plt.hist(recExcitonDisp, bins=10)
    plt.title('recExcitonDisp')
    plt.savefig('./recExcitonDisp.pdf')
    plt.clf()

    plt.hist(recExcitonTime, bins=10)
    plt.title('recExcitonTime')
    plt.savefig('./recExcitonTime.pdf')
    plt.clf()

    plt.hist(numberOfHops, bins=10)
    plt.title('numberOfHops')
    plt.savefig('./numberOfExcitonHops.pdf')
    plt.clf()

    print("XDE =", numberOfDissociations/parameterDict['minimumNumberOfPhotoinjections'])

    #plotConnections(excitonPath, chromophoreData.returnChromophoreList([1,8,6]), morphologyData.returnAAMorphology([1,8,6]))

    #print("Plotting the MSD of the exciton")
    #plt.figure()
    #MSD = list(np.array(excitonDisp)**2)
    #plt.plot(excitonTime, MSD, c='b')
    #fit = np.polyfit(excitonTime, MSD, 1)
    #print("Fit =", fit)
    #print("Diffusion Coefficient =", fit[0])
    #xFit = np.linspace(0, max(excitonTime), 100)
    #yFit = [(xVal*fit[0] + fit[1]) for xVal in xFit]
    #plt.plot(xFit, yFit, c='b')
    #plt.savefig('./excitonMSD.pdf')

    #print("Plotting the extraction data as a function of simulation time")
    #plt.figure()
    #plt.plot(convGlobalTime, convExtractions)
    #plt.savefig('./convergence.pdf')
    #return


if __name__ == "__main__":
    print("WOBBEY")
