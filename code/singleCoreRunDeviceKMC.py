import os
import sys
import copy
import random as R
import numpy as np
import heapq
import pickle
import subprocess as sp
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
epsilonNought = 8.85418782E-12     # m^{-3} kg^{-1} s^{4} A^{2}
ELECTRON = 0
HOLE = 1

# Global Variables
globalChromophoreData = []   # Contains all of the chromophoreLists for each morphology moiety for access by the hopping routines
globalMorphologyData = []    # Contains all of the AAMorphologyDicts for each morphology moiety
globalTime = 0               # The total simulation time that is required/updated all over the place
globalCarrierDict = {}       # A dictionary of all carriers in the system for the potential/recombination calculations
currentFieldValue = 0        # The field value calculated from the voltage this child process was given
numberOfExtractions = 0      # The total number of charges that have hopped out of the device through the `correct' contact


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
        self.prefactor = parameterDict['excitonPrefactor']
        self.numberOfHops = 0
        # NOTE For now, we'll just inject it randomly somewhere in the system
        self.initialDevicePosn = initialPosnDevice
        self.currentDevicePosn = copy.deepcopy(initialPosnDevice)
        self.initialChromophore = globalChromophoreData.returnRandomChromophore(initialPosnDevice)
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
                self.electronChromophore = globalChromophoreData.returnSpecificChromophore(self.currentDevicePosn, electronChromophoreID)
            elif self.currentChromophore.species == 'Acceptor':
                self.electronChromophore = self.currentChromophore
                # The hole chromophore is a randomly selected chromophore of the opposing type that is in range of the current one
                holeChromophoreID = R.choice(self.currentChromophore.dissociationNeighbours)[0]
                self.holeChromophore = globalChromophoreData.returnSpecificChromophore(self.currentDevicePosn, holeChromophoreID)
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
        # The prefactor included here is a bit of a bodge to try and get the mean-free paths of the excitons more in line with the 5nm of experiment. Possible citation: 10.3390/ijms131217019 (they don't do the simulation they just point out some limitations of FRET which assumes point-dipoles which doesn't necessarily work in all cases)
        if deltaEij <= 0:
            boltzmannFactor = 1
        else:
            boltzmannFactor = np.exp(-(elementaryCharge * deltaEij)/(kB * self.T))

        kFRET = self.prefactor * (1/self.lifetimeParameter) * (self.rF / rij)**6 * boltzmannFactor
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
        chromophoreList = globalChromophoreData.returnChromophoreList(self.currentDevicePosn)
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
            morphologyShape = [globalMorphologyData.returnAAMorphology(self.currentDevicePosn)[key] for key in ['lx', 'ly', 'lz']]
            remappedPosn = neighbourPosn + [neighbourRelativeImage[axis] * morphologyShape[axis] for axis in range(3)]
            rij = helperFunctions.calculateSeparation(currentChromoPosn, remappedPosn) * 1E-10
            #if neighbourRelativeImage != [0, 0, 0]:
            #    print(neighbourID, neighbourPosn, neighbourUnwrappedPosn, neighbourImage, neighbourRelativeImage, rij)
            # Note, separations are recorded in angstroems, so spin this down to metres.
            # Additionally, all of the energies are in eV currently, so convert them to J
            hopRate = self.calculateHopRate(rij, deltaEij * elementaryCharge)
            hopTime = self.determineHopTime(hopRate)
            # Keep track of the destination chromophore ID, the corresponding tau, and the relative image (to see if we're hopping over a boundary)
            hopTimes.append([globalChromophoreData.returnChromophoreList(self.currentDevicePosn)[self.currentChromophore.neighbours[neighbourIndex][0]], hopTime, neighbourRelativeImage])
        # Sort by ascending hop time
        hopTimes.sort(key = lambda x:x[1])
        # Only want to take the quickest hop if it is NOT going to hop outside the device (top or bottom - i.e. z axis is > shape[2] or < 0)
        # Other-axis periodic hops (i.e. x and y) are permitted, and will allow the exciton to loop round again.
        while len(hopTimes) > 0:
            # Get the hop destination cell
            newLocation = list(np.array(self.currentDevicePosn) + np.array(hopTimes[0][2]))
            # If the z axis component is > shape[2] or < 0, then forbid this hop by removing it from the hopTimes list.
            if (newLocation[2] >= globalChromophoreData.deviceArray.shape[2]) or (newLocation[2] < 0):
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
            newChromophore = globalChromophoreData.returnClosestChromophoreToPosition(self.currentDevicePosn, destinationPosition)
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
        # self.carrierType: Set carrier type to be == 0 if Electron and == 1 if Hole. This allows us to do quick arithmetic to get the signs correct in the potential calculations without having to burn through a ton of conditionals.
        if self.currentChromophore.species == 'Donor':
            self.carrierType = HOLE
            self.lambdaij = parameterDict['reorganisationEnergyDonor']
        elif self.currentChromophore.species == 'Acceptor':
            self.carrierType = ELECTRON
            self.lambdaij = parameterDict['reorganisationEnergyAcceptor']
        self.relativePermittivity = parameterDict['relativePermittivity']
        self.calculateBehaviour()

    def calculateBehaviour(self):
        [self.destinationChromophore, self.hopTime, self.destinationImage] = self.calculateHop()
        #try:
        #    [self.destinationChromophore, self.hopTime, self.destinationImage] = self.calculateHop()
        #except:
        #    [self.destinationChromophore, self.hopTime, self.destinationImage] = [None, None, None]

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
            destinationChromophore = globalChromophoreData.returnChromophoreList(self.currentDevicePosn)[self.currentChromophore.neighbours[neighbourIndex][0]]
            # Need to determine the relative image between the two (recorded in obtainChromophores) to check if the carrier is hopping out of this morphology cell and into an adjacent one
            neighbourRelativeImage = self.currentChromophore.neighbours[neighbourIndex][1]
            deltaEij = self.calculatePotential(destinationChromophore, neighbourRelativeImage, self.currentChromophore.neighboursDeltaE[neighbourIndex])
            # All of the energies are in eV currently, so convert them to J
            hopRate = self.calculateHopRate(self.lambdaij * elementaryCharge, transferIntegral * elementaryCharge, deltaEij * elementaryCharge)
            hopTime = self.determineHopTime(hopRate)
            # Keep track of the chromophoreID and the corresponding tau
            hopTimes.append([destinationChromophore, hopTime, neighbourRelativeImage])
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
            newChromophore = globalChromophoreData.returnClosestChromophoreToPosition(self.currentDevicePosn, destinationPosition)
            if (newChromophore == 'Top') or (newChromophore == 'Bottom'):
                # This carrier is hopping out of the active layer and into the contacts.
                # Firstly, ensure that it doesn't get queued up again
                [self.destinationChromophore, self.hopTime, self.destinationImage] = [None, None, None]
                self.removedTime = globalTime
                # Secondly, work out whether this is a `correct' hop (i.e. hole hopping to anode or electron hopping to cathode) that causes photovoltaic current.
                if newChromophore == 'Top':
                    # Leaving through top (anode)
                    if self.injectedFrom != 'Anode':
                        if self.carrierType == HOLE:
                            numberOfExtractions += 1
                        else:
                            numberOfExtractions -= 1
                    # Else (injected from anode), number of extractions doesn't change.
                else:
                    # Leaving through bottom (cathode)
                    if self.injectedFrom != 'Cathode':
                        if self.carrierType == ELECTRON:
                            numberOfExtractions += 1
                        else:
                            numberOfExtractions -= 1
                    # Else (injected from cathode), number of extractions doesn't change.
                print('EVENT: Charge Carrier Left Device! New number of extractions:', numberOfExtractions)
                input('Celebratory Pause...')
            else:
                self.currentChromophore = newChromophore

    def calculatePotential(self, destinationChromophore, neighbourRelativeImage, chromoEij):
        global globalCarrierDict

        # The potential will be calculated in eV
        currentAbsolutePosition = np.array(self.currentDevicePosn)*parameterDict['morphologyCellSize'] + (np.array(self.currentChromophore.posn) * 1E-10)
        destinationZ = (np.array(self.currentDevicePosn) + np.array(neighbourRelativeImage))*parameterDict['morphologyCellSize'] + (np.array(destinationChromophore.posn) * 1E-10)[2]
        zSep = destinationZ - currentAbsolutePosition[2]
        charge = (2 * self.carrierType) - 1
        totalPotential = (zSep * currentFieldValue * charge) + (chromoEij * charge)
        # Coulombic component!
        coulombicPotential = 0.0
        coulombConstant = 1.0 / (4 * np.pi * epsilonNought * self.relativePermittivity)
        for carrier in globalCarrierDict:
            # Each carrier feels coulombic effects from every other carrier in the system, along with the image charges induced in the top and bottom contacts, AND the image charges induced by the image charges.
            # The image charges are needed in order to keep the contacts as perfect metals (i.e. no internal electric field)
            # Check papers from Chris Groves and Ben Lyons for more details.
            carrierPosn = np.array(carrier.currentDevicePosn)*parameterDict['morphologyCellSize'] + (np.array(carrier.currentChromophore.posn) * 1E-10)
            if carrier.ID != self.ID:
                # Only consider the current carrier if it is not us!
                separation = np.linalg.norm(currentAboslutePosition - carrierPosn)
                coulombicPotential += ((2 * self.carrierType) - 1) * (coulombConstant / separation)
            # Now consider the image charges and the images of the images.
            # First do the top image charge
            deviceZSize = globalChromophoreData.deviceArray.shape[2] * parameterDict['morphologyCellSize']
            topImagePosn = copy.deepcopy(carrierPosn)
            topImagePosn[2] = deviceZSize + (deviceZSize - carrierPosn[2])
            separation = np.linalg.norm(currentAbsolutePosition - topImagePosn)
            coulombicPotential -= ((2 * self.carrierType) - 1) * (coulombConstant / separation)
            # Now do the bottom image charge
            bottomImagePosn = copy.deepcopy(carrierPosn)
            bottomImagePosn[2] = - carrierPosn[2]
            separation = np.linalg.norm(currentAbsolutePosition - bottomImagePosn)
            coulombicPotential -= ((2 * self.carrierType) - 1) * (coulombConstant / separation)
            # Now do the top image of the bottom image charge
            topImageOfBottomImagePosn = copy.deepcopy(carrierPosn)
            topImageOfBottomImagePosn[2] = (2 * deviceZSize) + carrierPosn[2]
            separation = np.linalg.norm(currentAbsolutePosition - topImageOfBottomImagePosn)
            coulombicPotential -= ((2 * self.carrierType) -1 ) * (coulombConstant / separation)
            # And finally the bottom image of the top image charge
            bottomImageOfTopImagePosn = copy.deepcopy(carrierPosn)
            bottomImageOfTopImagePosn[2] = - (deviceZSize + (deviceZSize - carrierPosn[2]))
            separation = np.linalg.norm(currentAbsolutePosition - bottomImageOfTopImagePosn)
            coulombicPotential -= ((2 * self.carrierType) - 1) * (coulombConstant / separation)
        print(coulombicPotential, zSep * currentFieldValue * charge, chromoEij * charge)
        raise SystemError("CHECK MAGNITUDES BEFORE CONTINUING")
        totalPotential += coulombicPotential
        return totalPotential


def plotHopDistance(distribution):
    plt.figure()
    plt.hist(distribution)
    plt.show()
    exit()


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


def execute(deviceArray, chromophoreData, morphologyData, parameterDict):
    # ---=== PROGRAMMER'S NOTE ===---
    # This `High Resolution' version of the code will permit chromophore-based hopping through the device, rather than just approximating the distribution. We have the proper resolution there, let's just use it!
    # ---=========================---
    global globalChromophoreData
    global globalMorphologyData
    global globalTime
    global currentFieldValue

    globalChromophoreData = chromophoreData
    globalMorphologyData = morphologyData

    # DEBUG
    #morphCellSizes = []
    #for moietyType in morphologyData.moietyDictionary.keys():
    #    morphCellSizes.append(morphologyData.moietyDictionary[moietyType].AAMorphologyDict['lx'])
    #print(np.average(morphCellSizes))
    #exit()

    # Need to initalise a bunch of variables
    eventQueue = []
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
                if injectedExciton.canDissociate is True:
                    print("EVENT: Exciton Dissociating!")
                    numberOfDissociations += 1
                    numberOfHops.append(injectedExciton.numberOfHops)
                    injectedElectron = carrier(carrierIndex, globalTime, injectedExciton.currentDevicePosn, injectedExciton.electronChromophore, 'Exciton', parameterDict)
                    globalCarrierDict[carrierIndex] = injectedElectron
                    carrierIndex += 1
                    injectedHole = carrier(carrierIndex, globalTime, injectedExciton.currentDevicePosn, injectedExciton.holeChromophore, 'Exciton', parameterDict)
                    # Add the hole to the carrier list for when we need to calc deltaE in the device
                    globalCarrierDict[carrierIndex] = injectedHole
                    carrierIndex += 1
                    # Now add both carriers' next hops to the KMC queue
                    heapq.heappush(eventQueue, (injectedElectron.hopTime, 'carrierHop', injectedElectron))
                    heapq.heappush(eventQueue, (injectedHole.hopTime, 'carrierHop', injectedHole))
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
            ## DEBUG CODE USED TO PLOT THE EXCITON ROUTE WITHIN A SINGLE CELL TO FIX THE PBCs
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
                    # Plop an electron down on the electronChromophore of the dissociating exciton
                    injectedElectron = carrier(carrierIndex, globalTime, hoppingExciton.currentDevicePosn, hoppingExciton.electronChromophore, 'Exciton', parameterDict)
                    # Add the electron to the carrier list for when we need to calc deltaE in the device
                    globalCarrierDict[carrierIndex] = injectedElectron
                    carrierIndex += 1
                    # Plop a hole down on the holeChromophore of the dissociating exciton
                    injectedHole = carrier(carrierIndex, globalTime, hoppingExciton.currentDevicePosn, hoppingExciton.holeChromophore, 'Exciton', parameterDict)
                    # Add the hole to the carrier list for when we need to calc deltaE in the device
                    globalCarrierDict[carrierIndex] = injectedHole
                    carrierIndex += 1
                    # Now add both carriers' next hops to the KMC queue
                    print(injectedElectron.hopTime)
                    print(injectedHole.hopTime)
                    raise SystemError("Check hoptimes are not none")
                    heapq.heappush(eventQueue, (injectedElectron.hopTime, 'carrierHop', injectedElectron))
                    heapq.heappush(eventQueue, (injectedHole.hopTime, 'carrierHop', injectedHole))
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
            hoppingCarrier = nextEvent[2]
            hoppingCarrier.performHop()
            if hoppingCarrier.removedTime is None:
                # As long as this carrier hasn't just been removed, recalculate its behaviour
                hoppingCarrier.calculateBehaviour()
            if hoppingCarrier.hopTime is None:
                # Carrier is either trapped with no eligible hops or has just been extracted
                if hoppingCarrier.removedTime is not None:
                    # Carrier has been extracted, so remove it from the carriers dictionary
                    globalCarrierDict.pop(hoppingCarrier.ID)
                # Else: Carrier is trapped, so we'll just leave it there so it affects the Coulombic landscape
            else:
                # Normal carrier hop, so requeue this carrier
                heapq.heappush(eventQueue, (hoppingCarrier.hopTime, 'carrierHop', hoppingCarrier))

        elif nextEvent[1] == 'recombine':
            raise SystemError("No recombine yet")

        else:
            print(eventQueue)
            raise SystemError("New Event is next in queue")

    ### THE FOLLOWING CODE WAS USED TO GENERATE THE ANALYSIS DATA USED IN THE MATTYSUMMARIES EXCITON DYNAMICS STUFF
    #print("Plotting Exciton Characteristics")
    #plt.figure()
    #plt.hist(dissExcitonDisp, bins=15)
    #plt.title('Displacement until Dissociation')
    #plt.xlabel('Displacement, nm')
    #plt.ylabel('Frequency, Arb. U.')
    #plt.savefig('./dissExcitonDisp.pdf')
    #plt.clf()

    #plt.hist(np.array(dissExcitonTime) * 1E9, bins=15)
    #plt.title('Time until Dissociation')
    #plt.xlabel('Time, ns')
    #plt.ylabel('Frequency, Arb. U.')
    #plt.savefig('./dissExcitonTime.pdf')
    #plt.clf()

    #plt.hist(recExcitonDisp, bins=15)
    #plt.title('Displacement until Recombination')
    #plt.xlabel('Displacement, nm')
    #plt.ylabel('Frequency, Arb. U.')
    #plt.savefig('./recExcitonDisp.pdf')
    #plt.clf()

    #plt.hist(np.array(recExcitonTime) * 1E10, bins=15)
    #plt.title('Time until Recombination')
    #plt.xlabel('Time, ns')
    #plt.ylabel('Frequency, Arb. U.')
    #plt.savefig('./recExcitonTime.pdf')
    #plt.clf()

    #plt.hist(numberOfHops, bins=15)
    #plt.title('numberOfHops')
    #plt.savefig('./numberOfExcitonHops.pdf')
    #plt.clf()

    #print("XDE =", numberOfDissociations/parameterDict['minimumNumberOfPhotoinjections'])
    #print("Mean/SD DissDisp =", np.mean(dissExcitonDisp), "+-", np.std(dissExcitonDisp)/np.sqrt(len(dissExcitonDisp)))
    #print("Mean/SD DissTime =", np.mean(dissExcitonTime), "+-", np.std(dissExcitonTime)/np.sqrt(len(dissExcitonTime)))
    #print("Mean/SD RecDisp =", np.mean(recExcitonDisp), "+-", np.std(recExcitonDisp)/np.sqrt(len(recExcitonDisp)))
    #print("Mean/SD RecTime =", np.mean(recExcitonTime), "+-", np.std(recExcitonTime)/np.sqrt(len(recExcitonTime)))

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
    KMCDirectory = sys.argv[1]
    CPURank = int(sys.argv[2])
    R.seed(int(sys.argv[3]))
    overwrite = False
    try:
        overwrite = bool(sys.argv[4])
    except:
        pass

    jobsFileName = KMCDirectory + '/KMCData_%02d.pickle' % (CPURank)
    deviceDataFileName = KMCDirectory.replace('/KMC', '/code/deviceData.pickle')

    with open(deviceDataFileName, 'rb') as pickleFile:
        [deviceArray, chromophoreData, morphologyData, parameterDict] = pickle.load(pickleFile)
    with open(jobsFileName, 'rb') as pickleFile:
        jobsToRun = pickle.load(pickleFile)
    logFile = KMCDirectory + '/KMClog_' + str(CPURank) + '.log'
    # Reset the log file
    with open(logFile, 'wb+') as logFileHandle:
        pass
    print('Found ' + str(len(jobsToRun)) + ' jobs to run.')
    #helperFunctions.writeToFile(logFile, ['Found ' + str(len(jobsToRun)) + ' jobs to run.'])
    # Set the affinities for this current process to make sure it's maximising available CPU usage
    currentPID = os.getpid()
    try:
        affinityJob = sp.Popen(['taskset', '-pc', str(CPURank), str(currentPID)], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
        # helperFunctions.writeToFile(logFile, affinityJob[0].split('\n')) #stdOut for affinity set
        # helperFunctions.writeToFile(logFile, affinityJob[1].split('\n')) #stdErr for affinity set
    except OSError:
        print("Taskset command not found, skipping setting of processor affinity...")
        #helperFunctions.writeToFile(logFile, ["Taskset command not found, skipping setting of processor affinity..."])

    # Begin the simulation
    execute(deviceArray, chromophoreData, morphologyData, parameterDict)
