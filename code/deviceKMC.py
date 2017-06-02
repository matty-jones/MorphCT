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
        try:
            deviceMoietyType = self.deviceArray[tuple(devicePosition)]
        except IndexError:
            # There is no device position that exists at these coordinates, and so we have a hop taking place that leaves the active layer
            # TODO THIS IS WHERE I GOT UP TO BEFORE EVAN

        for chromophore in moietyDictionary[deviceMoietyType].chromophoreList:
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
        # NOTE For now, we'll just inject it randomly somewhere in the system
        self.initialDevicePosn = initialPosnDevice
        self.currentDevicePosn = copy.deepcopy(initialPosnDevice)
        self.initialChromophore = chromophoreData.returnRandomChromophore(initialPosnDevice)
        self.currentChromophore = copy.deepcopy(self.currentChromophore)
        self.calculateBehaviour()

    def calculateBehaviour(self):
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
                holeChromophoID = R.choice(self.currentChromophore.dissociationNeighbours)[0]
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
        if globalTime > self.creationTime + self.recombinationTime:
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
            # Need to determine the relative image between the two (recorded in obtainChromophores). We will also need this to check if the exciton is hopping out of this morphology cell and into an adjacent one
            neighbourRelativeImage = self.currentChromophore.neighbours[neighbourIndex][1]
            # Morphology shape so we can work out the actual relative positions
            morphologyShape = [morphologyData.returnAAMorphology(self.currentDevicePosn)[key] for key in ['lx', 'ly', 'lz']]
            neighbourPosn += [neighbourRelativeImage[axis] * morphologyShape[axis] for axis in range(3)]
            rij = helperFunctions.calculateSeparation(currentChromoPosn, neighbourPosn) * 1E-10
            # Note, separations are recorded in angstroems, so spin this down to metres.
            # Additionally, all of the energies are in eV currently, so convert them to J
            hopRate = self.calculateHopRate(rij, deltaEij * elementaryCharge)
            hopTime = self.determineHopTime(hopRate)
            # Keep track of the destination chromophore ID, the corresponding tau, and the relative image (to see if we're hopping over a boundary)
            hopTimes.append([self.currentChromophore.neighbours[neighbourIndex][0], hopTime, neighbourRelativeImage])
        # Sort by ascending hop time
        hopTimes.sort(key = lambda x:x[1])
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
        deltaPosition = self.destinationPosition - initialPosition
        if self.destinationImage == [0, 0, 0]:
            # Exciton is not hopping over a boundary, so can simply update its current position
            self.currentChromophore = self.destinationChromophore
        else:
            # We're hopping over a boundary. Permit the hop with the already-calculated hopTime, but then find the closest chromophore to the destination position in the adjacent device cell.
            self.currentDevicePosn += self.destinationImage
            self.currentChromophore = chromophoreData.returnClosestChromophoreToPosition(self.currentDevicePosn, destinationPosition)
        # TODO: Now that the hop is complete, check to see if we can dissociate and then do so immediately if permitted
        self.canDissociate = self.checkDissociation()


class carrier:
    def __init__(self, index, globalTime, initialPosnDevice, initialChromophore, parameterDict):
        self.ID = carrierNo
        self.creationTime = globalTime
        self.removedTime = None
        self.initialDevicePosn = initialPosnDevice
        self.currentDevicePosn = copy.deepcopy(initialPosnDevice)
        self.initialChromophore = initialChromophore
        self.currentChromophore = copy.deepcopy(initialChromophore)
        self.T = parameterDict['systemTemperature']
        if self.currentChromophore.species == 'Donor':
            self.carrierType = 'Hole'
            self.lambdaij = parameterDict['reorganisationEnergyDonor']
        elif self.currentChromophore.species == 'Acceptor':
            self.carrierType = 'Electron'
            self.lambdaij = parameterDict['reorganisationEnergyAcceptor']
        self.calculateBehaviour()

    def calculateBehaviour(self):
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
            self.currentDevicePosn += self.destinationImage
            self.currentChromophore = chromophoreData.returnClosestChromophoreToPosition(self.currentDevicePosn, destinationPosn)
        # TODO: Incorporate the behaviour if the carrier is hopping out of the correct/incorrect electrode
        pass


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
    print("Decrementing time by", eventTime, "after event")
    globalTime += eventTime
    eventQueue = [(event[0] - eventTime, event[1], event[2]) for event in eventQueue]
    return eventQueue


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

    # Need to initalise a bunch of variables
    eventQueue = [None]
    # Need a carrier dictionary so that we can work out energetic landscape and remove the right carriers
    carriers = {}
    photoinjectionRate = calculatePhotoinjectionRate(parameterDict, deviceArray.shape)
    darkInjectionRate = calculateDarkInjectionRate(parameterDict, deviceArray.shape)
    outputCurrentDensity = []
    numberOfPhotoinjections = 0
    excitonIndex = 0
    carrierIndex = 0
    outputCurrentConverged = False
    print("\n\n ---=== MAIN KMC LOOP START ===---\n\n")
    # Main KMC loop
    while len(eventQueue) > 0:
        # Make sure that we run the loop at least once by appending None to eventQueue
        if eventQueue[0] is None:
            eventQueue.pop(0)
        # TODO Check whether the output current has converged, after we have had sufficient photoinjections
        if outputCurrentConverged is True:
            break

        # Consider photoinjections first
        photoinjectionTime = -np.log(R.random()) / photoinjectionRate
        heapq.heappush(eventQueue, (photoinjectionTime, 'photo', None))
        numberOfPhotoinjections += 1

        print("Event Queue Before Selection =", eventQueue)

        # Now find out what the next event is
        nextEvent = heapq.heappop(eventQueue)
        # Increment the global time and decrement all of the other times in the queue
        print("Event Queue after event has been popped", eventQueue)
        eventQueue = decrementTime(eventQueue, nextEvent[0])
        print("Event Queue after time has been decremented", eventQueue)

        # Execute the next behaviour (being sure to recalculate rates for any new particles that have appeared)
        if nextEvent[1] == 'photo':
            # Complete the event by injecting an exciton
            # First find an injection location. For now, this will just be somewhere random in the system
            randomDevicePosition = [R.randint(0, x-1) for x in deviceArray.shape]
            print("EVENT: Photoinjection into", randomDevicePosition, "which has type", deviceArray[tuple(randomDevicePosition)])
            injectedExciton = exciton(excitonIndex, globalTime, randomDevicePosition, parameterDict)
            # Push the exciton to the queue
            heapq.heappush(eventQueue, (injectedExciton.hopTime, 'excitonHop', injectedExciton))
            # Increment the exciton counter
            excitonIndex += 1
            print(eventQueue)
            input("Pause...")

        elif nextEvent[1] == 'dark':
            raise SystemError("No dark injection yet")

        elif nextEvent[1] == 'excitonHop':
            hoppingExciton = nextEvent[2]
            hoppingExciton.performHop()
            hoppingExciton.calculateBehaviour()
            # At this point, dissociate the exciton or remove it from the system if needed.
            if (hoppingExciton.canDissociate is True) or (hoppingExciton.hopTime is None):
                # Exciton needs to be removed. As we've already popped it from the queue, we just need to not queue it up again.
                if hoppingExciton.canDissociate is True:
                    injectedElectron = carrier(carrierIndex, globalTime, hoppingExciton.currentDevicePosn, hoppingExciton.electronChromophore, parameterDict)
                    # Add the electron to the carrier list for when we need to calc deltaE in the device
                    carriers[carrierIndex] = injectedElectron
                    carrierIndex += 1
                    injectedHole = carrier(carrierIndex, globalTime, hoppingExciton.currentDevicePosn, hoppingExciton.holeChromophore, parameterDict)
                    # Add the hole to the carrier list for when we need to calc deltaE in the device
                    carrier[carrierIndex] = injectedHole
                    carrierIndex += 1
                    # Now add both carriers' next hops to the KMC queue
                    heapq.heappush(eventQueue, (injectedElectron.hopTime, 'carrierHop', injectedElectron))
                    heapq.heappush(eventQueue, (injectedHole.hopTime, 'carrierHop', injectedHole))
                pass
            else:
                heapq.heappush(eventQueue, (hoppingExciton.hopTime, 'excitonHop', injectedExciton))

        elif nextEvent[1] == 'carrierHop':
            hoppingCarrier = nextEvent[2]
            hoppingCarrier.performHop()
            hoppingCarrier.calculateBehaviour()
            if hoppingCarrier.hopTime is not None:
                heapq.heappush(eventQueue, (hoppingCarrier.hopTime, 'carrierHop', hoppingCarrier))
            # Else: The carrier has become stuck (this should only be possible if we've crossed a device cell boundary and dropped the carrier on a chromophore with no corrections).
            # The best way to deal with this is to not queue it up, but leave it in the carrier dictionary so that it still affects the energetic landscape and can be recombined with.
            # Otherwise it will just sit here for the rest of the simulation and not do anything. Not sure if I like that.

        elif nextEvent[1] == 'recombine':
            raise SystemError("No recombine yet")

        else:
            print(eventQueue)
            raise SystemError("New Event is next in queue")

    return


if __name__ == "__main__":
    print("WOBBEY")
