import os
import sys
import copy
import random as R
import time as T
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
KMCIterations = 0
fastestEventAllowed = 1E-15
slowestEventAllowed = 1E-9


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
        self.prefactor = parameterDict['hoppingPrefactor']
        self.numberOfHops = 0
        # NOTE For now, we'll just inject it randomly somewhere in the system
        self.initialDevicePosn = initialPosnDevice
        self.currentDevicePosn = copy.deepcopy(initialPosnDevice)
        self.initialChromophore = globalChromophoreData.returnRandomChromophore(initialPosnDevice)
        self.currentChromophore = copy.deepcopy(self.initialChromophore)
        self.calculateBehaviour()
        if parameterDict['recordCarrierHistory'] is True:
            self.history = [[self.initialDevicePosn, self.initialChromophore.posn]]
        else:
            self.history = None

    def calculateBehaviour(self):
        # Set a flag to indicate whether exciton can dissociate or not
        self.canDissociate = self.checkDissociation()
        # Dissociate the exciton immediately after creation if it would be created at an interface
        if self.canDissociate is True:
            helperFunctions.writeToFile(logFile, ["self.canDissociate is True. self.currentChromophore.species = " +
                                                  str(self.currentChromophore.species)])
            # Determine all potential dissociation options, randomly select one, plop a hole on the donor and electron on the acceptor, then remove this exciton from the system by updating its removedTime
            if self.currentChromophore.species == 'Donor':
                holeChromophore = globalChromophoreData.returnSpecificChromophore(self.currentDevicePosn, self.currentChromophore.ID)
                # The electron chromophore is a randomly selected chromophore of the opposing type that is in range of the current one
                electronChromophoreID = R.choice(self.currentChromophore.dissociationNeighbours)[0]
                electronChromophore = globalChromophoreData.returnSpecificChromophore(self.currentDevicePosn, electronChromophoreID)
                if (self.currentDevicePosn not in holeChromophore.occupied) and (self.currentDevicePosn not in electronChromophore.occupied):
                    self.holeChromophore = holeChromophore
                    self.electronChromophore = electronChromophore
                else:
                    helperFunctions.writeToFile(logFile, ["Debug: Cannot dissociate after all. currentDevicePosn = " +
                                                          repr(self.currentDevicePosn) + " holeChromophore.occupied = " +
                                                          str(holeChromophore.occupied) + "electronChromophore.occupied = " +
                                                          str(electronChromophore.occupied)])
                    self.canDissociate = False
            elif self.currentChromophore.species == 'Acceptor':
                electronChromophore = globalChromophoreData.returnSpecificChromophore(self.currentDevicePosn, self.currentChromophore.ID)
                # The hole chromophore is a randomly selected chromophore of the opposing type that is in range of the current one
                holeChromophoreID = R.choice(self.currentChromophore.dissociationNeighbours)[0]
                holeChromophore = globalChromophoreData.returnSpecificChromophore(self.currentDevicePosn, holeChromophoreID)
                if (self.currentDevicePosn not in holeChromophore.occupied) and (self.currentDevicePosn not in electronChromophore.occupied):
                    self.holeChromophore = holeChromophore
                    self.electronChromophore = electronChromophore
                else:
                    helperFunctions.writeToFile(logFile, ["Debug: Cannot dissociate after all. currentDevicePosn = " +
                                                          repr(self.currentDevicePosn) + " holeChromophore.occupied = " +
                                                          str(holeChromophore.occupied) + "electronChromophore.occupied = " +
                                                          str(electronChromophore.occupied)])
                    self.canDissociate = False
            if self.canDissociate is True:
                # Notify execute() that this exciton should not be queued up again by setting self.hopTime == None
                [self.destinationChromophore, self.hopTime, self.destinationImage] = [None, None, None]
                self.removeTime = globalTime
                return
        # If we're not instantaneously dissociating, check the recombination time
        if globalTime >= self.creationTime + self.recombinationTime:
            # Exciton has run out of time and is now recombining, so set its next hopTime to None
            [self.destinationChromophore, self.hopTime, self.destinationImage] = [None, None, None]
            self.removedTime = globalTime
            return
        # Otherwise, calculate the fastest hop from the current chromophore
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
            rij = helperFunctions.calculateSeparation(np.array(currentChromoPosn), np.array(remappedPosn)) * 1E-10
            #if neighbourRelativeImage != [0, 0, 0]:
            #    print(neighbourID, neighbourPosn, neighbourUnwrappedPosn, neighbourImage, neighbourRelativeImage, rij)
            # Note, separations are recorded in angstroems, so spin this down to metres.
            # Additionally, all of the energies are in eV currently, so convert them to J
            hopRate = calculateFRETHopRate(self.prefactor, self.lifetimeParameter, self.rF, rij, deltaEij * elementaryCharge, self.T)
            hopTime = determineEventTau(hopRate, 'exciton-hop')
            # Keep track of the destination chromophore ID, the corresponding tau, and the relative image (to see if we're hopping over a boundary)
            if hopTime is not None:
                hopTimes.append([globalChromophoreData.returnSpecificChromophore(self.currentDevicePosn, self.currentChromophore.neighbours[neighbourIndex][0]), hopTime, neighbourRelativeImage])
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
            targetDevicePosn = list(np.array(self.currentDevicePosn) + np.array(self.destinationImage))
            newChromophore = globalChromophoreData.returnClosestChromophoreToPosition(targetDevicePosn, destinationPosition)
            if (newChromophore == 'Top') or (newChromophore == 'Bottom'):
                # This exciton is hopping out of the active layer and into the contacts.
                # Ensure it doesn't get queued up again
                [self.destinationChromophore, self.hopTime, self.destinationImage] = [None, None, None]
                self.removedTime = globalTime
                # TODO This carrier has left the simulation volume so we now need to ensure we remove it from the carrier dictionary so it's not included in any of the Coulombic calculations
            elif (newChromophore == 'Out of Bounds'):
                # Ignore this hop, do a different one instead
                helperFunctions.writeToFile(logFile, ["Out of bounds Exciton hop discarded"])
            else:
                self.currentDevicePosn = list(np.array(self.currentDevicePosn) + np.array(self.destinationImage))
                self.currentChromophore = newChromophore
        self.numberOfHops += 1
        self.canDissociate = self.checkDissociation()
        if self.history is not None:
            self.history.append([self.currentDevicePosn, self.currentChromophore.posn])


class carrier:
    def __init__(self, index, globalTime, initialPosnDevice, initialChromophore, injectedFrom, parameterDict, injectedOntoSite = None):
        self.ID = index
        self.creationTime = globalTime
        self.removedTime = None
        self.initialDevicePosn = initialPosnDevice
        self.currentDevicePosn = initialPosnDevice
        self.initialChromophore = initialChromophore
        self.currentChromophore = initialChromophore
        self.injectedFrom = injectedFrom
        self.injectedOntoSite = injectedOntoSite
        self.prefactor = parameterDict['hoppingPrefactor']
        self.recombining = False
        self.recombiningWith = None
        self.T = parameterDict['systemTemperature']
        self.wrapxy = parameterDict['wrapDeviceXY']
        # self.carrierType: Set carrier type to be == 0 if Electron and == 1 if Hole. This allows us to do quick arithmetic to get the signs correct in the potential calculations without having to burn through a ton of conditionals.
        if self.currentChromophore.species == 'Donor':
            self.carrierType = HOLE
            self.lambdaij = parameterDict['reorganisationEnergyDonor']
        elif self.currentChromophore.species == 'Acceptor':
            self.carrierType = ELECTRON
            self.lambdaij = parameterDict['reorganisationEnergyAcceptor']
        self.relativePermittivity = parameterDict['relativePermittivity']
        if parameterDict['recordCarrierHistory'] is True:
            self.history = [[self.initialDevicePosn, initialChromophore.posn]]
        else:
            self.history = None

    def calculateBehaviour(self):
        try:
            [self.destinationChromophore, self.hopTime, self.relativeImage, self.destinationImage] = self.calculateHop()
            # Update the destination chromophore (if it's a chromophore) so that it's marked as occupied
            if not isinstance(self.destinationChromophore, str):
                globalChromophoreData.returnSpecificChromophore(self.destinationImage, self.destinationChromophore.ID).occupied.append(self.destinationImage)
        except ValueError:
            [self.destinationChromophore, self.hopTime, self.relativeImage, self.destinationImage] = [None, None, None, None]

    def calculateHop(self):
        # Determine the hop times to all possible neighbours
        hopTimes = []
        # Obtain the reorganisation energy in J (from eV in the parameter file)
        for neighbourIndex, transferIntegral in enumerate(self.currentChromophore.neighboursTI):
            # Ignore any hops with a NoneType transfer integral (usually due to an ORCA error), or zero
            if (transferIntegral is None) or (transferIntegral < 1E-10):
                continue
            neighbourChromophore = globalChromophoreData.returnSpecificChromophore(self.currentDevicePosn, self.currentChromophore.neighbours[neighbourIndex][0])
            # The destination chromophore will be the actual chromophore we end up on (i.e. not neighbour if we hop across a boundary)
            destinationChromophore = neighbourChromophore
            # Need to determine the relative image between the two (recorded in obtainChromophores) to check if the carrier is hopping out of this morphology cell and into an adjacent one
            neighbourRelativeImage = self.currentChromophore.neighbours[neighbourIndex][1]
            destinationImage = list(np.array(self.currentDevicePosn) + np.array(neighbourRelativeImage))
            # If we're hopping out of this cell and into a new one, make sure that it's in-range if we're not going to wrap it
            skipThisNeighbour = False
            if self.wrapxy is False:
                # Skip this neighbour if we're going to hop out of the device along the non-electrode axes
                if neighbourRelativeImage != [0, 0, 0]:
                    for axisNo, val in enumerate(destinationImage[:-1]):
                        if (val >= globalChromophoreData.deviceArray.shape[axisNo]) or (val < 0):
                            skipThisNeighbour = True
                            break
                    # Don't need to perform this check if we've already satisfied the first skip condition
                    if skipThisNeighbour is False:
                        destinationChromophore = globalChromophoreData.returnClosestChromophoreToPosition(destinationImage, neighbourChromophore.posn)
                        if destinationChromophore == "Out of Bounds":
                            continue
                    else:
                        continue
            # Make sure the destination chromophore is of the correct type and is not occupied, if it is a chromophore and not a string (saying 'Top' or 'Bottom' if it's leaving the device)
            if not isinstance(destinationChromophore, str):
                if (destinationImage in destinationChromophore.occupied) or (destinationChromophore.species != self.currentChromophore.species):
                    continue
            # Otherwise, we're good to go. Calculate the hop as previously (but a hop to the neighbourChromophore)
            deltaEij = self.calculateDeltaE(neighbourChromophore, neighbourRelativeImage, self.currentChromophore.neighboursDeltaE[neighbourIndex])
            # All of the energies (EXCEPT EIJ WHICH IS ALREADY IN J) are in eV currently, so convert them to J
            hopRate = calculateMarcusHopRate(self.prefactor, self.lambdaij * elementaryCharge, transferIntegral * elementaryCharge, deltaEij, self.T)
            #print(hopRate, self.prefactor, self.lambdaij * elementaryCharge, transferIntegral * elementaryCharge, deltaEij)
            #print(self.currentChromophore.ID, "to", destinationChromophore.ID)
            hopTime = determineEventTau(hopRate, 'carrier-hop')
            if hopTime is not None:
                # Keep track of the chromophoreID and the corresponding tau
                hopTimes.append([destinationChromophore, hopTime, neighbourRelativeImage, destinationImage, self.prefactor, self.lambdaij * elementaryCharge, transferIntegral * elementaryCharge, deltaEij])
                if hopTime > 1E20:
                    helperFunctions.writeToFile(logFile, ["---=== WARNING, EXTREMELY LONG CARRIER HOP ===---", repr(hopTimes[-1])])
        # Sort by ascending hop time
        hopTimes.sort(key = lambda x:x[1])
        # Take the quickest hop
        if len(hopTimes) > 0:
            if hopTimes[0][1] > 1E98:
                helperFunctions.writeToFile(logFile, [repr(hopTimes[0]), repr(globalCarrierDict)])
                for ID, carrier in globalCarrierDict.items():
                    helperFunctions.writeToFile(logFile, [str(ID), repr(carrier.currentDevicePosn),
                                                          str(carrier.currentChromophore.ID)])
            return hopTimes[0][:4]
        else:
            # We are trapped here, so just return in order to raise the calculateBehaviour exception
            return []

    def performHop(self):
        global numberOfExtractions
        global KMCIterations
        global globalTime
        # Unset the current chromophore's occupation status as we're about to hop away.
        # The destination should have already been updated
        globalChromophoreData.returnSpecificChromophore(self.currentDevicePosn, self.currentChromophore.ID).occupied.remove(self.currentDevicePosn)
        if (self.destinationChromophore == 'Top') or (self.destinationChromophore == 'Bottom'):
                # This carrier is hopping out of the active layer and into the contacts.
                # Firstly, work out whether this is a `correct' hop (i.e. hole hopping to anode or electron hopping to cathode) that causes photovoltaic current.
                if self.destinationChromophore == 'Top':
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
                if self.carrierType == ELECTRON:
                    helperFunctions.writeToFile(logFile, ["EVENT: Electron Left out of " + str(self.destinationChromophore) +
                                                          " of Device! New number of extractions: " + str(numberOfExtractions) +
                                                          " after " + str(KMCIterations) + " iterations (globalTime = " +
                                                          str(globalTime) + ")"])
                else:
                    helperFunctions.writeToFile(logFile, ["EVENT: Hole Left out of " + str(self.destinationChromophore) +
                                                          " of Device! New number of extractions: " + str(numberOfExtractions) +
                                                          " after " + str(KMCIterations) + " iterations (globalTime = " +
                                                          str(globalTime) + ")"])
                # Now ensure that it doesn't get queued up again
                [self.destinationChromophore, self.hopTime, self.relativeImage, self.destinationImage] = [None, None, None, None]
                self.removedTime = globalTime
        else:
            self.currentDevicePosn = self.destinationImage
            self.currentChromophore = self.destinationChromophore
        if self.history is not None:
            self.history.append([self.currentDevicePosn, self.currentChromophore.posn])

    def calculateDeltaE(self, destinationChromophore, neighbourRelativeImage, chromoEij):
        # DeltaEij has 3 main components: 1) the energetic disorder (difference in HOMO/LUMO levels), 
        # 2) the field within the device, and 3) the Coulombic effect from nearby charges
        deltaEij = 0.0  # Report this in J
        # 1) Energetic Disorder
        deltaEij += (chromoEij * elementaryCharge)
        # 2) Field within the device
        currentAbsolutePosition = np.array(self.currentDevicePosn) * parameterDict['morphologyCellSize'] + (np.array(self.currentChromophore.posn) * 1E-10)
        destinationAbsolutePosition = ((np.array(self.currentDevicePosn) + np.array(neighbourRelativeImage)) * parameterDict['morphologyCellSize'] + (np.array(destinationChromophore.posn) * 1E-10))
        # Field has negative sign because device is flipped with anode at +Z and Cathode at 0
        zSep = - (destinationAbsolutePosition[2] - currentAbsolutePosition[2])
        charge = elementaryCharge * ((2 * self.carrierType) - 1)
        deltaEij += (zSep * currentFieldValue * charge)
        # 3) Difference in Coulombic Potential at dest compared to origin
        originCoulomb, recombineFlag, recombineID = calculateCoulomb(currentAbsolutePosition, self.ID, self.carrierType, carrierIsRecombining=self.recombining)
        destinationCoulomb, dummy1, dummy2 = calculateCoulomb(destinationAbsolutePosition, self.ID, self.carrierType)
        deltaEij += (destinationCoulomb - originCoulomb)
        if (recombineFlag is True) and (self.recombining is False):
            # Need to update the recombining details
            self.recombining = recombineFlag
            self.recombiningWith = recombineID
        return deltaEij


class injectSite:
    def __init__(self, devicePosn, chromophore, injectRate, electrode):
        self.devicePosn = devicePosn
        self.chromophore = chromophore
        self.injectRate = injectRate
        self.electrode = electrode
        self.calculateInjectTime()

    def calculateInjectTime(self):
        #self.injectTime = np.float64(1.0 / self.injectRate)
        self.injectTime = determineEventTau(self.injectRate, self.electrode + '-injection')


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


def decrementTime(eventQueue, eventTime):
    global globalTime
    # A function that increments the global time by whatever the time this event is, and
    # then decrements all of the remaining times in the queue.
    globalTime += eventTime
    eventQueue = [(event[0] - eventTime, event[1], event[2]) for event in eventQueue]
    return eventQueue


def plotConnections(excitonPath, chromophoreList, AADict, outputDir):
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
    fileName = outputDir + '3d_exciton.pdf'
    plt.show()
    plt.savefig(fileName)
    print("Figure saved as " + fileName)


def calculateCoulomb(absolutePosition, selfID, selfCarrierType, carrierIsRecombining=None):
    global globalCarrierDict
    coulombicPotential = 0.0
    coulombConstant = 1.0 / (4 * np.pi * epsilonNought * parameterDict['relativePermittivity'])
    canRecombineHere = False
    recombiningWith = None
    for carrierID, carrier in globalCarrierDict.items():
        # Each carrier feels coulombic effects from every other carrier in the system, along with the image charges induced in the top and bottom contacts, AND the image charges induced by the image charges.
        # The image charges are needed in order to keep a constant field at the interface, correcting for the non-periodic Neumann-like boundary conditions
        # Check papers from Chris Groves (Marsh, Groves and Greenham2007 and beyond), Ben Lyons (2011/2012), and van der Holst (2011) for more details.
        carrierPosn = np.array(carrier.currentDevicePosn)*parameterDict['morphologyCellSize'] + (np.array(carrier.currentChromophore.posn) * 1E-10)
        if carrierID != selfID:
            # Only consider the current carrier if it is not us!
            separation = helperFunctions.calculateSeparation(absolutePosition, carrierPosn)
            #print("Separation between", selfID, "and", carrier.ID, "=", separation)
            coulombicPotential += coulombConstant * ((elementaryCharge * ((2 * selfCarrierType) - 1)) * (elementaryCharge * ((2 * carrier.carrierType) - 1)) / separation)
            # I'm also going to use this opportunity (as we're iterating over all carriers in the system) to see if we're close enough to any to recombine.
            if carrierIsRecombining is False:
                # Only do this if we're not currently recombining with something
                if (separation <= parameterDict['coulombCaptureRadius']) and (selfCarrierType != carrier.carrierType):
                    # If carriers are within 1nm of each other, then assume that they are within the Coulomb capture radius and are about to recombine
                    helperFunctions.writeToFile(logFile, ["CARRIER " + str(selfID) + " CLOSE ENOUGH " +
                                                          "TO RECOMBINE WITH " + str(carrierID)])
                    canRecombineHere = True
                    recombiningWith = carrierID
                    # We don't necessarily need to update the carrier.recombining, since all of the behaviour will be determined in the main program by checking hoppingCarrier.recombining, but just in case, we can update the carrier too as well as self.
                    globalCarrierDict[carrierID].recombining = True
                    globalCarrierDict[carrierID].recombiningWith = selfID
        # TODO: NOTE: DEBUG: TURNING OFF IMAGE CHARGES TO SEE IF I CAN WORK OUT WHY EVERYTHING IS FAILING
        ## Now consider the image charges and the images of the images.
        ## Direct images have the opposing charge to the actual carrier, so the potential contains (+q * -q), or -(q**2)
        ## Images of the images have the opposing charge to the image of the carrier, which is the same charge as the actual carrier
        ## Therefore, images of images are either (+q * +q) or (-q * -q) = +(q**2)
        ## First do the top image charge
        #deviceZSize = globalChromophoreData.deviceArray.shape[2] * parameterDict['morphologyCellSize']
        #topImagePosn = copy.deepcopy(carrierPosn)
        #topImagePosn[2] = deviceZSize + (deviceZSize - carrierPosn[2])
        #separation = helperFunctions.calculateSeparation(absolutePosition, topImagePosn)
        #print("Separation between", selfID, "and top image charge =", separation)
        #coulombicPotential -= (elementaryCharge**2) * (coulombConstant / separation)
        ## Now do the bottom image charge
        #bottomImagePosn = copy.deepcopy(carrierPosn)
        #bottomImagePosn[2] = - carrierPosn[2]
        #separation = helperFunctions.calculateSeparation(absolutePosition, bottomImagePosn)
        #print("Separation between", selfID, "and bottom image charge =", separation)
        #coulombicPotential -= (elementaryCharge**2) * (coulombConstant / separation)
        ## Now do the top image of the bottom image charge
        #topImageOfBottomImagePosn = copy.deepcopy(carrierPosn)
        #topImageOfBottomImagePosn[2] = (2 * deviceZSize) + carrierPosn[2]
        #separation = helperFunctions.calculateSeparation(absolutePosition, topImageOfBottomImagePosn)
        #print("Separation between", selfID, "and top image of bottom image charge =", separation)
        #coulombicPotential += (elementaryCharge**2) * (coulombConstant / separation)
        ## And finally the bottom image of the top image charge
        #bottomImageOfTopImagePosn = copy.deepcopy(carrierPosn)
        #bottomImageOfTopImagePosn[2] = - (deviceZSize + (deviceZSize - carrierPosn[2]))
        #separation = helperFunctions.calculateSeparation(absolutePosition, bottomImageOfTopImagePosn)
        #print("Separation between", selfID, "and bottom image of top image charge =", separation, "\n")
        #coulombicPotential += (elementaryCharge**2) * (coulombConstant / separation)
    return coulombicPotential, canRecombineHere, recombiningWith


def calculateDarkCurrentInjections(deviceArray, parameterDict):
    # NOTE Had a nightmare with Marcus hopping here, and I'm not convinced it makes sense to have hops from metals
    # be described wrt reorganisation energies and transfer integrals. Instead, I'm going to use a more generic
    # MA hopping methodology here. This therefore needs a prefactor, which is defined in the parameter file.
    # This is another variable to calibrate, but for now I've set it such that we get a reasonable hopping rate
    # for dark-injection hops (somewhere in the middle of our normal hopping-rate distribution (i.e. ~1E12))
    # Get the device shape so that we know which cells to inject into
    deviceShape = deviceArray.shape
    # Calculate the important energy levels
    bandgap = parameterDict['acceptorLUMO'] - parameterDict['donorHOMO']
    electronInjectBarrier = parameterDict['acceptorLUMO'] - parameterDict['cathodeWorkFunction']
    holeInjectBarrier = parameterDict['anodeWorkFunction'] - parameterDict['donorHOMO']
    validCathodeInjSites = 0
    validAnodeInjSites = 0
    cathodeInjectRatesData = []
    anodeInjectRatesData = []
    cathodeInjectWaitTimes = []
    anodeInjectWaitTimes = []
    # Consider electron-injecting electrode (cathode) at bottom of device first:
    zVal = 0
    for xVal in range(deviceShape[0]):
        for yVal in range(deviceShape[1]):
            cathodeInjectChromophores = []
            AAMorphology = globalMorphologyData.returnAAMorphology([xVal, yVal, zVal])
            morphologyChromophores = globalChromophoreData.returnChromophoreList([xVal, yVal, zVal])
            for chromophore in morphologyChromophores:
                # Find all chromophores that are between 5 and 10 Ang from the bottom 10 Ang of the device cell
                if (chromophore.posn[2] <= -(AAMorphology['lz'] / 2.0) + 10) and (chromophore.posn[2] >= -(AAMorphology['lz'] / 2.0) + 5) and (len([_ for _ in chromophore.neighboursTI if (_ is not None) and (_ > 1E-5)]) > 0):
                #if (chromophore.posn[2] <= -(AAMorphology['lz'] / 2.0) + 10) and (sum(1 for _ in filter(None.__ne__, chromophore.neighboursDeltaE)) > 0):
                    cathodeInjectChromophores.append([chromophore, 1E-10 * (chromophore.posn[2] - (-(AAMorphology['lz'] / 2.0)))])
                    validCathodeInjSites += 1
            for [chromophore, separation] in cathodeInjectChromophores:
                if chromophore.species == 'Acceptor':
                    # Injecting an electron from the cathode (easy)
                    #print("EASY INJECT", electronInjectBarrier, "+", chromophore.LUMO - parameterDict['acceptorLUMO'])
                    #print("Chromo LUMO =", chromophore.LUMO)
                    #print("Target LUMO =", parameterDict['acceptorLUMO'])
                    deltaE = elementaryCharge * (electronInjectBarrier)
                    #deltaE = elementaryCharge * (electronInjectBarrier - (chromophore.LUMO - parameterDict['acceptorLUMO']))
                    #print("Easy ELECTRON Inject from Cathode =", deltaE / elementaryCharge)
                    #input("PAUSE...")
                else:
                    # Injecting a hole from the cathode (hard)
                    #print("Difficult INJECT", bandgap, "-", electronInjectBarrier, "-", (chromophore.LUMO - parameterDict['acceptorLUMO']))
                    #print("Chromo LUMO =", chromophore.LUMO)
                    #print("Target LUMO =", parameterDict['donorHOMO'])
                    # Here, we make the assumption that the energy difference between the chromophore and literature HOMOs
                    # would be the same as the energy difference between the chromophore and literature LUMOs for injection
                    # (this is because we don't record the Donor LUMO or the Acceptor HOMO)
                    deltaE = elementaryCharge * (bandgap - electronInjectBarrier)
                    #deltaE = elementaryCharge * (bandgap - electronInjectBarrier - (chromophore.HOMO - parameterDict['donorHOMO']))
                    #print("Difficult HOLE inject from Cathode =", deltaE / elementaryCharge)
                    #input("PAUSE...")
                #print("CATHODE DELTA E", deltaE, chromophore.species)
                injectRate = calculateMillerAbrahamsHopRate(parameterDict['MAPrefactor'], separation, parameterDict['MALocalisationRadius'], deltaE, parameterDict['systemTemperature'])
                #print("CATHODE", injectRate)
                cathodeInjectRatesData.append(injectRate)
                # Create inject site object
                site = injectSite([xVal, yVal, zVal], chromophore, injectRate, 'Cathode')
                injectTime = site.injectTime
                cathodeInjectWaitTimes = pushToQueue(cathodeInjectWaitTimes, (injectTime, 'cathode-injection', site))
    # Now consider hole-injecting electrode at top of device (anode):
    #input("PAUSE...")
    zVal = deviceShape[2]-1
    for xVal in range(deviceShape[0]):
        for yVal in range(deviceShape[1]):
            anodeInjectChromophores = []
            AAMorphology = globalMorphologyData.returnAAMorphology([xVal, yVal, zVal])
            morphologyChromophores = globalChromophoreData.returnChromophoreList([xVal, yVal, zVal])
            estimatedTransferIntegral = estimateTransferIntegral(morphologyChromophores) * elementaryCharge
            for chromophore in morphologyChromophores:
                # Find all chromophores that are between 5 and 10 Ang from the top of the device cell
                if (chromophore.posn[2] >= (AAMorphology['lz'] / 2.0) - 10) and (chromophore.posn[2] <= (AAMorphology['lz'] / 2.0) - 5) and (len([_ for _ in chromophore.neighboursTI if (_ is not None) and (_ > 1E-5)]) > 0):
                #if (chromophore.posn[2] >= (AAMorphology['lz'] / 2.0) - 10) and (sum(1 for _ in filter(None.__ne__, chromophore.neighboursDeltaE)) > 0):
                    anodeInjectChromophores.append([chromophore, 1E-10 * ((AAMorphology['lz'] / 2.0) - chromophore.posn[2])])
                    validAnodeInjSites += 1
            for [chromophore, separation] in anodeInjectChromophores:
                if chromophore.species == 'Acceptor':
                    # Injecting an electron from the anode (hard)
                    deltaE = elementaryCharge * (bandgap - holeInjectBarrier)
                    #deltaE = elementaryCharge * (bandgap - holeInjectBarrier - (chromophore.LUMO - parameterDict['acceptorLUMO']))
                    #print("Hard ELECTRON Inject from Anode =", deltaE / elementaryCharge)
                else:
                    # Injecting a hole from the anode (easy)
                    deltaE = elementaryCharge * (holeInjectBarrier)
                    #deltaE = elementaryCharge * (holeInjectBarrier + (chromophore.HOMO - parameterDict['donorHOMO']))
                    #print("Easy HOLE Inject from Anode =", deltaE / elementaryCharge)
                #print("ANODE DELTA E", deltaE, chromophore.species)
                injectRate = calculateMillerAbrahamsHopRate(parameterDict['MAPrefactor'], separation, parameterDict['MALocalisationRadius'], deltaE, parameterDict['systemTemperature'])
                anodeInjectRatesData.append(injectRate)
                # Create inject site object
                site = injectSite([xVal, yVal, zVal], chromophore, injectRate, 'Anode')
                injectTime = site.injectTime
                anodeInjectWaitTimes = pushToQueue(anodeInjectWaitTimes, (injectTime, 'anode-injection', site))
    #plt.figure()
    #plt.hist(injectRates, bins = np.logspace(9, 12, 30))
    #plt.gca().set_xscale('log')
    #print("Valide Cathode Injection Sites =", validCathodeInjSites)
    #print("Valide Anode Injection Sites =", validAnodeInjSites)
    #print("Mean =", np.mean(injectRates), "STD =", np.std(injectRates))
    #print("Max =", max(injectRates), "Min =", min(injectRates))
    #plt.show()
    #exit()
    cathodeInjectRate = np.mean(cathodeInjectRatesData)
    anodeInjectRate = np.mean(anodeInjectRatesData)
    return cathodeInjectRate, anodeInjectRate, cathodeInjectWaitTimes, anodeInjectWaitTimes


def getNextDarkEvent(queue, electrode):
    # A function that pops the next event from a queue (cathode or anode injection queues), decrements the time, requeues the injectSite with a new time, pushes to the queue and returns the next event and the updated queue
    nextEvent = heapq.heappop(queue)
    queue = decrementTime(queue, nextEvent[0])
    requeuedEvent = copy.deepcopy(nextEvent[2])
    requeuedEvent.calculateInjectTime()
    queue = pushToQueue(queue, (requeuedEvent.injectTime, electrode + '-injection', requeuedEvent))
    return nextEvent, queue


def estimateTransferIntegral(chromophoreList):
    # Iterate over all non-periodic chromophore hopping events in the system
    TIs = []
    for chromophore in chromophoreList:
        for index, chromo in enumerate(chromophore.neighbours):
            if chromo[1] != [0, 0, 0]:
                continue
            TI = chromophore.neighboursTI[index]
            if TI is not None:
                TIs.append(TI)
    return np.average(TIs)


    ## NOTE The following code doesn't work very well so I scrapped it and used a different technique to get the TI:
    ## Fits a straight line to a particular chromophore's neighbour transfer integrals as a function of separation, then interpolates this trend for a given startingZPosn to estimate the TI given a particular separation.
    ## This is important for calculating an estimated TI for dark-current injected carriers from an electrode to a particular chromophore, which can be a variable-range hop to any chromophore within 1nm of the electrode
    #separations = []
    #transferIntegrals = []
    #currentPosn = copy.deepcopy(chromophore.posn)
    #currentPosn[2] = startingZPosn
    #for index, chromo in enumerate(chromophore.neighbours):
    #    if chromo[1] != [0, 0, 0]:
    #        # Ignore hops that are not in this simulation volume
    #        continue
    #    chromoID = chromo[0]
    #    separation = helperFunctions.calculateSeparation(currentPosn, chromophoreList[chromoID].posn)
    #    separations.append(separation)
    #    transferIntegrals.append(chromophore.neighboursTI[index])
    #if len(separations) > 1:
    #    fit = np.polyfit(separations, transferIntegrals, 1)
    #else:
    #    # Only one datapoint so set the fit to be y = 0x + TI
    #    fit = [0, transferIntegrals[0]]
    #interpolatedTI = (fit[0] * startingZPosn) + fit[1]
    #if interpolatedTI < 0:
    #    raise SystemError("INTERPOLATION FAIL")
    #input("PAUSE")
    #return interpolatedTI


def calculateMillerAbrahamsHopRate(prefactor, separation, radius, deltaEij, T):
    kij = prefactor * np.exp(-2 * separation/radius)
    if deltaEij > 0:
        kij *= np.exp(-deltaEij / (kB * T))
    #print("kij =", kij)
    #print("prefactor =", prefactor)
    #print("separation =", separation)
    #print("radius =", radius)
    #input("PAUSE")
    return kij


def calculateMarcusHopRate(prefactor, lambdaij, Tij, deltaEij, T):
    ## Semiclassical Marcus Hopping Rate Equation
    #print("\nLambdaij =", lambdaij)
    #print("Tij =", Tij)
    #print("DeltaEij =", deltaEij)
    #print("Marcus components:")
    #print(((2 * np.pi) / hbar) * (Tij ** 2))
    #print(np.sqrt(1.0 / (4 * lambdaij * np.pi * kB * T)))
    #print(np.exp(-((deltaEij + lambdaij)**2) / (4 * lambdaij * kB * T)))
    #print((deltaEij + lambdaij) ** 2)
    #print(4 * lambdaij * kB * T)
    kij = prefactor * ((2 * np.pi) / hbar) * (Tij ** 2) * np.sqrt(1.0 / (4 * lambdaij * np.pi * kB * T)) * np.exp(-((deltaEij + lambdaij)**2) / (4 * lambdaij * kB * T))
    #print(kij)
    #input("PAUSE...")
    return kij


def calculateFRETHopRate(prefactor, lifetimeParameter, rF, rij, deltaEij, T):
    # Foerster Transport Hopping Rate Equation
    # The prefactor included here is a bit of a bodge to try and get the mean-free paths of the excitons more in line with the 5nm of experiment. Possible citation: 10.3390/ijms131217019 (they do not do the simulation they just point out some limitations of FRET which assumes point-dipoles which does not necessarily work in all cases)
    if deltaEij <= 0:
        boltzmannFactor = 1
    else:
        boltzmannFactor = np.exp(-(elementaryCharge * deltaEij)/(kB * T))

    kFRET = prefactor * (1/lifetimeParameter) * (rF / rij)**6 * boltzmannFactor
    return kFRET


def determineEventTau(rate, eventType):
    # Use the KMC algorithm to determine the wait time to this event
    if rate != 0:
        counter = 0
        while True:
            if counter == 100:
                if 'hop' in eventType:
                    return None
                else:
                    if 'injection' not in eventType:
                        helperFunctions.writeToFile(logFile, ["Attempted 100 times to obtain a '" + str(eventType) +
                                                              "'-type event timescale within the tolerances: " +
                                                              str(fastestEventAllowed) + " <= tau < " +
                                                              str(slowestEventAllowed) + " with the given rate " +
                                                              str(rate) + " all without success.",
                                                              "Permitting the event anyway with the next random number."])
            x = R.random()
            # Ensure that we do not get exactly 0.0 or 1.0, which would break our logarithm
            if (x == 0) or (x == 1):
                continue
            tau = - np.log(x) / rate
            if ((tau > fastestEventAllowed) and (tau < slowestEventAllowed)) or (counter == 100):
                break
            else:
                counter += 1
    else:
        # If rate == 0, then make the event wait time extremely long
        tau = np.float64(1E99)
    return tau


def pushToQueue(queue, event):
    event = list(event)
    try:
        event[0] = np.float64(event[0])
    except:
        helperFunctions.writeToFile(logFile, ["Tried to push an event to the queue that does not have a" +
                                              "np.float64 as the first element (and so can't be queued)",
                                              "Queue = " + repr(queue), "Event = " + repr(event),
                                              "Terminating..."])
        raise KeyboardInterrupt
    if (event[0] == np.float64(1E99)) and (event[1] != 'photo'):
        helperFunctions.writeToFile(logFile, ["---=== TRIED TO QUEUE EVENT WITH CRAZY LONG WAIT TIME ===---",
                                              "Event = " + str(event), "This carrier has completed " +
                                              str(len(event[2].history)) + " hops."])
        try:
            event[2].__dict__.pop('history')
        except KeyError:
            pass
        helperFunctions.writeToFile(logFile, ["Instance = " + repr(event[2].__dict__)])
        try:
            helperFunctions.writeToFile(logFile, ["Current Chromophore = " + repr(event[2].currentChromophore.__dict__),
                                                  "Destination Chromophore = " + repr(event[2].destinationChromophore.__dict__)])
        except AttributeError:
            pass
        helperFunctions.writeToFile(logFile, ["Terminating..."])
        for index, carrier in globalCarrierDict.items():
            print(carrier.currentDevicePosn, carrier.currentChromophore.posn)
        for carrier1ID, carrier1 in globalCarrierDict.items():
            for carrier2ID, carrier2 in globalCarrierDict.items():
                if carrier1ID >= carrier2ID:
                    continue
                carrier1posn = (78 * np.array(carrier1.currentDevicePosn)) + np.array(carrier1.currentChromophore.posn)
                carrier2posn = (78 * np.array(carrier2.currentDevicePosn)) + np.array(carrier2.currentChromophore.posn)
                print(carrier1ID, carrier2ID, helperFunctions.calculateSeparation(carrier1posn, carrier2posn))
        raise KeyboardInterrupt
    event = tuple(event)
    heapq.heappush(queue, event)
    return queue


def plotEventTimeDistribution(eventLog, outputDir, fastest, slowest):
    fig = plt.figure()
    plt.hist([eventLog[_] for _ in eventLog.keys()], bins=np.logspace(int(np.floor(np.log10(fastest))), int(np.ceil(np.log10(slowest))), 10), color=['r', 'g', 'c', 'm', 'b', 'y'], label=eventLog.keys(), linewidth=0)
    plt.legend(loc = 1, prop = {'size':6})
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.xlabel(r'$\mathrm{\tau}$ (s)')
    plt.ylabel('Freq (Arb. U.)')
    plt.xticks([1E-16, 1E-14, 1E-12, 1E-10, 1E-8, 1E-6, 1E-5, 1E-4])
    fileName = outputDir + 'EventTimeDist.pdf'
    plt.savefig(fileName)
    print('Event time distribution saved as ' + fileName)


def gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


def plotCarrierTrajectories(allCarriers, parameterDict, deviceArray, outputDir):
    combinationsToPlot = {'Cathode':[], 'Anode':[]}
    popList = []
    for index, carrier in enumerate(allCarriers):
        # First add all the excitons
        if 'rF' in carrier.__dict__:
            # We know this is an exciton
            combinationsToPlot[carrier.ID] = [carrier]
            popList.append(index)
    for index in sorted(popList, reverse=True):
        allCarriers.pop(index)
    for carrier in allCarriers:
        # Only electrons and holes remaining
        try:
            combinationsToPlot[carrier.injectedFrom].append(carrier)
        except KeyError:
            pass
    for injectSource, carriers in combinationsToPlot.items():
        plot3DTrajectory(injectSource, carriers, parameterDict, deviceArray, outputDir)


def plot3DTrajectory(injectSource, carriersToPlot, parameterDict, deviceArray, outputDir):
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    [xLen, yLen, zLen] = 1E10 * (np.array(deviceArray.shape) * parameterDict['morphologyCellSize'])
    # The conversion is needed to get the box in ang
    # Draw boxlines
    # Varying X
    ax.plot([0, xLen], [0, 0], [0, 0], c = 'k', linewidth = 1.0)
    ax.plot([0, xLen], [yLen, yLen], [0, 0], c = 'k', linewidth = 1.0)
    ax.plot([0, xLen], [0, 0], [zLen, zLen], c = 'k', linewidth = 1.0)
    ax.plot([0, xLen], [yLen, yLen], [zLen, zLen], c = 'k', linewidth = 1.0)
    # Varying Y
    ax.plot([0, 0], [0, yLen], [0, 0], c = 'k', linewidth = 1.0)
    ax.plot([xLen, xLen], [0, yLen], [0, 0], c = 'k', linewidth = 1.0)
    ax.plot([0, 0], [0, yLen], [zLen, zLen], c = 'k', linewidth = 1.0)
    ax.plot([xLen, xLen], [0, yLen], [zLen, zLen], c = 'k', linewidth = 1.0)
    # Varying Z
    ax.plot([0, 0], [0, 0], [0, zLen], c = 'k', linewidth = 1.0)
    ax.plot([0, 0], [yLen, yLen], [0, zLen], c = 'k', linewidth = 1.0)
    ax.plot([xLen, xLen], [0, 0], [0, zLen], c = 'k', linewidth = 1.0)
    ax.plot([xLen, xLen], [yLen, yLen], [0, zLen], c = 'k', linewidth = 1.0)

    if (injectSource is 'Anode') or (injectSource is 'Cathode'):
        carrierString = injectSource
    else:
        # Exciton
        carrierString = 'Exciton_%04d' % (carriersToPlot[0].ID)
    for carrier in carriersToPlot:
        if 'rF' in carrier.__dict__:
            color = 'g'
        elif carrier.carrierType == ELECTRON:
            color = 'b'
        elif carrier.carrierType== HOLE:
            color = 'r'
        for hopIndex, hop in enumerate(carrier.history[:-1]):
            # Note that conversion factors are needed as hop[0] * morphCellSize is in m, and hop[1] is in ang.
            # Additionally, we need to add half of the morphCellSize to whatever hop[0] is as the origin is in the centre of the box 
            currentPosn = (1E10 * ((np.array(hop[0]) + np.array([0.5, 0.5, 0.5])) * parameterDict['morphologyCellSize'])) + np.array(hop[1])
            nextHop = carrier.history[hopIndex + 1]
            nextPosn = (1E10 * ((np.array(nextHop[0]) + np.array([0.5, 0.5, 0.5])) * parameterDict['morphologyCellSize'])) + np.array(nextHop[1])
            ax.plot([currentPosn[0], nextPosn[0]], [currentPosn[1], nextPosn[1]], [currentPosn[2], nextPosn[2]], c = color, linewidth = 0.5)
    fileName = outputDir + carrierString + '_traj.pdf'
    plt.savefig(fileName)
    print("Figure saved as " + fileName)
    plt.close()


def execute(deviceArray, chromophoreData, morphologyData, parameterDict, voltageVal, timeLimit):
    # ---=== PROGRAMMER'S NOTE ===---
    # This `High Resolution' version of the code will permit chromophore-based hopping through the device, rather than just approximating the distribution. We have the proper resolution there, let's just use it!
    # ---=========================---
    global globalChromophoreData
    global globalMorphologyData
    global globalTime
    global currentFieldValue
    global KMCIterations
    global fastestEventAllowed
    global slowestEventAllowed
    global logFile

    globalChromophoreData = chromophoreData
    globalMorphologyData = morphologyData
    if parameterDict['fastestEventAllowed'] is not None:
        fastestEventAllowed = parameterDict['fastestEventAllowed']
    if parameterDict['slowestEventAllowed'] is not None:
        slowestEventAllowed = parameterDict['slowestEventAllowed']
    # Given a voltage, the field value corresponding to it = ((bandgap - el_inj_barrier - ho_inj_barrier) / z-extent) 
    currentFieldValue = (
        # Bandgap:
        ((parameterDict['acceptorLUMO'] - parameterDict['donorHOMO']) -
        # Electron Inject Barrier:
        (parameterDict['acceptorLUMO'] - parameterDict['cathodeWorkFunction']) -
        # Hole Inject Barrier:
        (parameterDict['anodeWorkFunction'] - parameterDict['donorHOMO'])) /
        # Z-extent:
        (deviceArray.shape[2] * parameterDict['morphologyCellSize']))
    outputFiguresDir = parameterDict['outputDeviceDir'] + '/' + parameterDict['deviceMorphology'] + '/figures/' + str(voltageVal) +'/'

    # DEBUG
    slowestEvent = 0
    fastestEvent = 2E99
    eventTimes = []
    #morphCellSizes = []
    #for moietyType in morphologyData.moietyDictionary.keys():
    #    morphCellSizes.append(morphologyData.moietyDictionary[moietyType].AAMorphologyDict['lx'])
    #print(np.average(morphCellSizes))
    #exit()

    # Need to initalise a bunch of variables
    eventQueue = []
    photoinjectionRate = calculatePhotoinjectionRate(parameterDict, deviceArray.shape)
    outputCurrentDensity = []
    numberOfPhotoinjections = 0
    numberOfCathodeInjections = 0
    numberOfAnodeInjections = 0
    excitonIndex = 0
    carrierIndex = 0
    outputCurrentConverged = False
    recombiningCarrierIDs = []

    if parameterDict['recordCarrierHistory'] is True:
        allCarriers = []



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
    eventLog = {'photo':[], 'cathode-injection':[], 'anode-injection':[], 'excitonHop':[], 'carrierHop':[], 'recombine':[]}
    alreadyPrinted = False
    outputPrintStatement = False


    # As the morphology is not changing, we can calculate the dark inject rates and locations at the beginning and then not have to do it again, just re-queue up the injectSites as we use them by calling site(calculateInjectTime)
    cathodeInjectRate, anodeInjectRate, cathodeInjectQueue, anodeInjectQueue = calculateDarkCurrentInjections(deviceArray, parameterDict)

    #darkInjectRateDist = [_[0] for _ in darkInjectRatesData]
    #print("Dark injection rate stats:")
    #print("Mean =", np.mean(darkInjectRateDist), "STD =", np.std(darkInjectRateDist), "Min =", min(darkInjectRateDist), "Max =", max(darkInjectRateDist))
    #plt.figure()
    #plt.hist(darkInjectRateDist, bins = np.logspace(2, 10, 20))
    #plt.gca().set_xscale('log')
    #plt.show()
    #exit()

    print("\n\n ---=== MAIN KMC LOOP START ===---\n\n")
    t0 = T.time()
    # Main KMC loop
    # Put this in a try-except to plot on keyboardInterrupt
    try:
        while True:
            KMCIterations += 1
            t1 = T.time()
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
                if (numberOfPhotoinjections > parameterDict['minimumNumberOfPhotoinjections']):
                    # If we've gone over the minimum, then perform the convergence check to see if the output current has converged
                    # TODO Check convergence, but for now, let's just see what the curve looks like
                    break

            # Break if current has converged
            if outputCurrentConverged is True:
                break

            # Break if less than an hour before SLURM will kill the job
            if (timeLimit is not None) and (t1 > t0 + timeLimit - 3600):
                break

            if len(eventQueue) == 0:
                # Either the simulation has just started, or everything just recombined/left the device
                # Therefore, queue up one photoinjection and one dark current injection to get us rolling.
                # These are the only available kinetic starting points for the simulation. Everything else stems from these.
                # Photoinjection
                photoinjectionTime = determineEventTau(photoinjectionRate, 'photoinjection')
                eventQueue = pushToQueue(eventQueue, (photoinjectionTime, 'photo', None))
                numberOfPhotoinjections += 1
                # Dark Injection:
                # We now split up dark injection to be separate between the anode and the cathode, to ensure detail balance.
                # We will use the cathodeInjectRate (the mean of the distribution) and the anodeInjectRate to determine when to perform the injection
                # Then, we will pop the first event of the cathodeInjectWaitTimes and anodeInjectWaitTimes as required, executing it IMMEDIATELY and decrementing the time in the rest of the queue.
                # Then, we can re-queue the inject site.
                # This means that the wait times in cathodeInjectWaitTimes and anodeInjectWaitTimes are actually only used to determine the order of the chromophores onto which we inject from each electrode, and do not increment the globalTime.
                # The global time is incremented according to cathodeInjectionTime and anodeInjectionTime instead.
                # First, deal with the cathode injection
                cathodeInjectionTime = determineEventTau(cathodeInjectRate, 'cathode-injection')
                # Sort out the cathodeInjectQueue by getting the next inject site, and requeueing it with a new time to update the queue.
                nextCathodeEvent, cathodeInjectQueue = getNextDarkEvent(cathodeInjectQueue, 'Cathode')
                # Now push the cathode event to the main queue
                eventQueue = pushToQueue(eventQueue, (cathodeInjectionTime, 'cathode-injection', nextCathodeEvent[2]))

                # Now, deal with the anode injection
                anodeInjectionTime = determineEventTau(anodeInjectRate, 'anode-injection')
                # Sort out the anodeInjectQueue by popping the next inject site, decrementing the queue and re-queueing the inject site
                nextAnodeEvent, anodeInjectQueue = getNextDarkEvent(anodeInjectQueue, 'Anode')
                # Now push the anode event to the main queue
                eventQueue = pushToQueue(eventQueue, (anodeInjectionTime, 'anode-injection', nextAnodeEvent[2]))

            if int(t1 - t0)%10 == 0:
                outputPrintStatement = True
            else:
                outputPrintStatement = False
                alreadyPrinted = False
            if outputPrintStatement and not alreadyPrinted:
                helperFunctions.writeToFile(logFile, ["Current runtime = " + str(int(t1-t0)) + "s, with " +
                                                      str(len(eventQueue)) + " events currently in the queue and " +
                                                      str(len(globalCarrierDict.keys())) +
                                                      " carriers currently in the system. Currently completed " +
                                                      str(KMCIterations) + " iterations and simulated " +
                                                      str(globalTime) + "s."])
                alreadyPrinted = True


            # Now find out what the next event is
            nextEvent = heapq.heappop(eventQueue)
            # Increment the global time and decrement all of the other times in the queue
            #print("Event Queue after event has been popped", eventQueue)
            eventQueue = decrementTime(eventQueue, nextEvent[0])
            #print("Event Queue after time has been decremented", eventQueue)

            # DEBUG
            if nextEvent[0] > slowestEvent:
                slowestEvent = nextEvent[0]
            if nextEvent[0] < fastestEvent:
                fastestEvent = nextEvent[0]
            eventTimes.append(nextEvent[0])
            eventLog[nextEvent[1]].append(nextEvent[0])


            # Execute the next behaviour (being sure to recalculate rates for any new particles that have appeared)
            if nextEvent[1] == 'photo':
                # Complete the event by injecting an exciton
                # First find an injection location. For now, this will just be somewhere random in the system
                randomDevicePosition = [R.randint(0, x-1) for x in deviceArray.shape]
                helperFunctions.writeToFile(logFile, ["EVENT: Photoinjection #" + str(numberOfPhotoinjections) +
                                                      " into " + repr(randomDevicePosition) + " (which has type " +
                                                      repr(deviceArray[tuple(randomDevicePosition)]) + ") after " +
                                                      str(KMCIterations) + " iterations (globalTime = " +
                                                      str(globalTime) + ")"])
                injectedExciton = exciton(excitonIndex, globalTime, randomDevicePosition, parameterDict)
                if (injectedExciton.canDissociate is True) or (injectedExciton.hopTime is None):
                    # Injected onto either a dissociation site or a trap site
                    if injectedExciton.canDissociate is True:
                        helperFunctions.writeToFile(logFile, ["EVENT: Exciton Dissociating immediately"])
                        numberOfDissociations += 1
                        numberOfHops.append(injectedExciton.numberOfHops)
                        # Create the carrier instances, but don't yet calculate their behaviour (need to add them to the carrier list before we can calculate the energetics)
                        # Also add the carriers to the carrier dictionary for when we need to calc deltaE in the device
                        injectedElectron = carrier(carrierIndex, globalTime, injectedExciton.currentDevicePosn, injectedExciton.electronChromophore, injectedExciton.ID, parameterDict)
                        globalCarrierDict[carrierIndex] = injectedElectron
                        carrierIndex += 1
                        injectedHole = carrier(carrierIndex, globalTime, injectedExciton.currentDevicePosn, injectedExciton.holeChromophore, injectedExciton.ID, parameterDict)
                        globalCarrierDict[carrierIndex] = injectedHole
                        carrierIndex += 1
                        # Update the injected chromophores to be marked as occupied
                        globalChromophoreData.returnSpecificChromophore(injectedExciton.currentDevicePosn, injectedExciton.electronChromophore.ID).occupied.append(injectedExciton.currentDevicePosn)
                        globalChromophoreData.returnSpecificChromophore(injectedExciton.currentDevicePosn, injectedExciton.holeChromophore.ID).occupied.append(injectedExciton.currentDevicePosn)
                        # Now determine the behaviour of both carriers, and add their next hops to the KMC queue
                        injectedElectron.calculateBehaviour()
                        if injectedElectron.hopTime is not None:
                            eventQueue = pushToQueue(eventQueue, (injectedElectron.hopTime, 'carrierHop', injectedElectron))
                        injectedHole.calculateBehaviour()
                        if injectedHole.hopTime is not None:
                            eventQueue = pushToQueue(eventQueue, (injectedHole.hopTime, 'carrierHop', injectedHole))
                    else:
                        # Injected onto a site with no connections, so this exciton will eventually die
                        helperFunctions.writeToFile(logFile, ["EVENT: Exciton Recombining immediately"])
                        numberOfRecombinations += 1
                        numberOfHops.append(injectedExciton.numberOfHops)
                else:
                    # Hopping permitted, push the exciton to the queue.
                    # Push the exciton to the queue
                    eventQueue = pushToQueue(eventQueue, (injectedExciton.hopTime, 'excitonHop', injectedExciton))
                    # Increment the exciton counter
                    excitonIndex += 1
                # A photoinjection has just occured, so now queue up a new one
                photoinjectionTime = determineEventTau(photoinjectionRate, 'photoinjection')
                eventQueue = pushToQueue(eventQueue, (photoinjectionTime, 'photo', None))
                numberOfPhotoinjections += 1

            elif (nextEvent[1] == 'cathode-injection') or (nextEvent[1] == 'anode-injection'):
                injectSite = nextEvent[2]
                if injectSite.electrode == 'Cathode':
                    numberOfInjections = numberOfCathodeInjections
                else:
                    numberOfInjections = numberOfAnodeInjections
                injectChromophore = globalChromophoreData.returnSpecificChromophore(injectSite.devicePosn, injectSite.chromophore.ID)
                if (injectSite.devicePosn not in injectChromophore.occupied):
                    helperFunctions.writeToFile(logFile, ["EVENT: Dark Current injection from the " +
                                                          str(injectSite.electrode) + " #" +
                                                          str(numberOfInjections) + " into " +
                                                          repr(injectSite.devicePosn) + " (which has type " +
                                                          repr(deviceArray[tuple(injectSite.devicePosn)]) +
                                                          ") chromophore number " + str(injectSite.chromophore.ID) +
                                                          " after " + str(KMCIterations) + " iterations (globalTime = " +
                                                          str(globalTime) + ")"])
                    # Inject the carrier
                    injectedCarrier = carrier(carrierIndex, globalTime, injectSite.devicePosn, injectSite.chromophore, injectSite.electrode, parameterDict, injectedOntoSite = injectSite)
                    globalCarrierDict[carrierIndex] = injectedCarrier
                    carrierIndex += 1
                    # Update the chromophore occupation
                    helperFunctions.writeToFile(logFile, ["PREVIOUS OCCUPATION = " + repr(injectChromophore.occupied)])
                    injectChromophore.occupied.append(injectSite.devicePosn)
                    helperFunctions.writeToFile(logFile, ["UPDATED ON INJECTION = " + repr(injectChromophore.occupied)])
                    # Determine the injected carrier's next hop and queue it
                    injectedCarrier.calculateBehaviour()
                    if (injectedCarrier.hopTime is not None) and (injectedCarrier.hopTime > 1):
                        helperFunctions.writeToFile(logFile, ["DARK INJECTION LED TO ELECTRON WITH CRAZY HOPTIME"])
                        for carrierFromList in globalCarrierDict.values():
                            print(carrierFromList.currentDevicePosn, carrierFromList.currentChromophore.posn)
                        helperFunctions.writeToFile(logFile, [str(injectedCarrier.currentChromophore.ID)])
                        helperFunctions.writeToFile(logFile, [repr(injectedCarrier.__dict__)])
                        exit()
                    if injectedCarrier.hopTime is not None:
                        eventQueue = pushToQueue(eventQueue, (injectedCarrier.hopTime, 'carrierHop', injectedCarrier))
                else:
                    helperFunctions.writeToFile(logFile, ["CANNOT INJECT HERE, DESTINATION CHROMOPHORE IS ALREADY OCCUPIED"])

                # Now determine the next DC event and queue it
                if injectSite.electrode == 'Cathode':
                    nextCathodeEvent, cathodeInjectQueue = getNextDarkEvent(cathodeInjectQueue, 'Cathode')
                    eventQueue = pushToQueue(eventQueue, (cathodeInjectionTime, 'cathode-injection', nextCathodeEvent[2]))
                    numberOfCathodeInjections += 1
                if injectSite.electrode == 'Anode':
                    nextAnodeEvent, anodeInjectQueue = getNextDarkEvent(anodeInjectQueue, 'Anode')
                    eventQueue = pushToQueue(eventQueue, (anodeInjectionTime, 'anode-injection', nextAnodeEvent[2]))
                    numberOfAnodeInjections += 1

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
                    if parameterDict['recordCarrierHistory'] is True:
                        allCarriers.append(hoppingExciton)

                    if hoppingExciton.canDissociate is True:
                        helperFunctions.writeToFile(logFile, ["EVENT: Exciton Dissociating after " +
                                                              str(KMCIterations) + " iterations (globalTime = " +
                                                              str(globalTime) + ")"])
                        numberOfDissociations += 1
                        numberOfHops.append(injectedExciton.numberOfHops)
                        # Create the carrier instances, but don't yet calculate their behaviour (need to add them to the carrier list before we can calculate the energetics)
                        # Also add the carriers to the carrier dictionary for when we need to calc deltaE in the device
                        # Plop an electron down on the electronChromophore of the dissociating exciton, and a hole on the holeChromophore
                        injectedElectron = carrier(carrierIndex, globalTime, hoppingExciton.currentDevicePosn, hoppingExciton.electronChromophore, hoppingExciton.ID, parameterDict)
                        globalCarrierDict[carrierIndex] = injectedElectron
                        carrierIndex += 1
                        injectedHole = carrier(carrierIndex, globalTime, hoppingExciton.currentDevicePosn, hoppingExciton.holeChromophore, hoppingExciton.ID, parameterDict)
                        globalCarrierDict[carrierIndex] = injectedHole
                        carrierIndex += 1
                        # Update the injected chromophores to be marked as occupied
                        globalChromophoreData.returnSpecificChromophore(hoppingExciton.currentDevicePosn, hoppingExciton.electronChromophore.ID).occupied.append(hoppingExciton.currentDevicePosn)
                        globalChromophoreData.returnSpecificChromophore(hoppingExciton.currentDevicePosn, hoppingExciton.holeChromophore.ID).occupied.append(hoppingExciton.currentDevicePosn)
                        # Add to the allCarriers list for plotting
                        if parameterDict['recordCarrierHistory'] is True:
                            allCarriers += [injectedElectron, injectedHole]
                        # Now determine the behaviour of both carriers, and add their next hops to the KMC queue
                        injectedElectron.calculateBehaviour()
                        if injectedElectron.hopTime is not None:
                            eventQueue = pushToQueue(eventQueue, (injectedElectron.hopTime, 'carrierHop', injectedElectron))
                        injectedHole.calculateBehaviour()
                        if injectedHole.hopTime is not None:
                            eventQueue = pushToQueue(eventQueue, (injectedHole.hopTime, 'carrierHop', injectedHole))
                    else:
                        helperFunctions.writeToFile(logFile, ["EVENT: Exciton Recombining after " + str(KMCIterations) +
                                                              " iterations (globalTime = " + str(globalTime) + ")"])
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
                    eventQueue = pushToQueue(eventQueue, (hoppingExciton.hopTime, 'excitonHop', injectedExciton))

            elif nextEvent[1] == 'carrierHop':
                hoppingCarrier = nextEvent[2]
                # Check that the carrier is still in the carrier dictionary. If it's not, then it has already been removed from the system and we can safely ignore it
                if hoppingCarrier.ID not in globalCarrierDict.keys():
                    continue
                #print(hoppingCarrier.carrierType, hoppingCarrier.currentChromophore.species, hoppingCarrier.destinationChromophore.species, hoppingCarrier.currentDevicePosn, "type:", globalMorphologyData.returnDeviceMoietyType(hoppingCarrier.currentDevicePosn))
                #input("Check this is right. CarrierType = 0 is electron, 1 is hole")
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
                    if parameterDict['recordCarrierHistory'] is True:
                        allCarriers.append(hoppingCarrier)
                else:
                    # Normal carrier hop, so requeue this carrier
                    eventQueue = pushToQueue(eventQueue, (hoppingCarrier.hopTime, 'carrierHop', hoppingCarrier))
                # Check if we're eligible to recombine
                if (hoppingCarrier.recombining is True) and (hoppingCarrier.ID not in recombiningCarrierIDs) and (hoppingCarrier.recombiningWith not in recombiningCarrierIDs):
                    recombiningCarrierIDs.append(hoppingCarrier.ID)
                    recombiningCarrierIDs.append(hoppingCarrier.recombiningWith)
                    recombinationTime = determineEventTau(parameterDict['recombinationRate'], 'carrier-recombination')
                    eventQueue = pushToQueue(eventQueue, (recombinationTime, 'recombine', hoppingCarrier))
                ## If this carrier was injected, requeue this injection site in the darkInjectQueue to allow another injection
                ## Only re-queueing up the dark inject site after the carrier has hopped prevents double injection onto the same
                ## chromophore, by waiting until the carrier hops away first.
                ## TODO, what if the carrier hops to a chromophore that was about to be injected into? We don't check for this, so
                ## could still get double occupancy!
                #if hoppingCarrier.injectedOntoSite is not None:
                #    newDarkInject = copy.deepcopy(injectSite)
                #    newDarkInject.calculateInjectTime()
                #    heapq.heappush(darkInjectQueue, (newDarkInject.injectTime, 'dark', newDarkInject))

            elif nextEvent[1] == 'recombine':
                helperFunctions.writeToFile(logFile, ["EVENT: Carrier Recombination Check after " + str(KMCIterations) +
                                                      " iterations (globalTime = " + str(globalTime) + ")"])
                # A recombination event is about to occur. At this point, we should check if the carrier and its recombination partner are still in range.
                try:
                    carrier1 = nextEvent[2]
                    carrier1Posn = np.array(carrier1.currentDevicePosn) * parameterDict['morphologyCellSize'] + (np.array(carrier1.currentChromophore.posn) * 1E-10)
                    carrier2 = globalCarrierDict[carrier1.recombiningWith]
                    carrier2Posn = np.array(carrier2.currentDevicePosn) * parameterDict['morphologyCellSize'] + (np.array(carrier2.currentChromophore.posn) * 1E-10)
                    separation = helperFunctions.calculateSeparation(carrier2Posn, carrier1Posn)
                    recombiningCarrierIDs.remove(carrier2.ID)
                except KeyError:
                    # The second carrier is missing from the simulation (already extracted), so set the separation to be large
                    separation = 1E99
                recombiningCarrierIDs.remove(carrier1.ID)
                if separation <= parameterDict['coulombCaptureRadius']:
                    helperFunctions.writeToFile(logFile, [str(separation) + " <= " + str(parameterDict['coulombCaptureRadius'])])
                    helperFunctions.writeToFile(logFile, ["EVENT: Carrier Recombination Succeeded after " + str(KMCIterations) +
                                                          " iterations (globalTime = " + str(globalTime) + ")"])
                    # Carriers are in range, so recombine them
                    carrier1.removedTime = globalTime
                    carrier2.removedTime = globalTime
                    numberOfRecombinations += 1
                    globalCarrierDict.pop(carrier1.ID)
                    globalCarrierDict.pop(carrier2.ID)
                    if parameterDict['recordCarrierHistory'] is True:
                        allCarriers += [carrier1, carrier2]
                else:
                    helperFunctions.writeToFile(logFile, [str(separation) + " >" + str(parameterDict['coulombCaptureRadius'])])
                    helperFunctions.writeToFile(logFile, ["EVENT: Carrier Recombination Failed after " + str(KMCIterations) +
                                                          " iterations (globalTime = " + str(globalTime) + ")"])
                    # Carriers are no longer in range, so the recombination fails. Update their recombination flags
                    carrier1.recombining = False
                    carrier1.recombiningWith = None
                    if separation != 1E99:
                        # If recombining carrier was already extracted, then skip this
                        carrier2.recombining = False
                        carrier2.recombiningWith = None

            else:
                print(eventQueue)
                raise SystemError("New Event is next in queue")

    except KeyboardInterrupt:
        time = T.time() - t0
        print("KILL COMMAND RECIEVED...")
        print("Plotting output graphs before terminating...")
        helperFunctions.writeToFile(logFile, ["---=== RESULTS FROM CPU RANK " + sys.argv[2] + " ===---",
                                              "Run terminated after " + str(KMCIterations) +
                                              " iterations (globalTime = " + str(globalTime) + ") after " +
                                              str(time) + " seconds",
                                              "Number of Photoinjections = " + str(numberOfPhotoinjections),
                                              "Number of Cathode Injections = " + str(numberOfCathodeInjections),
                                              "Number of Anode Injections = " + str(numberOfAnodeInjections),
                                              "Number of Dissociations = " + str(numberOfDissociations),
                                              "Number of Recombinations = " + str(numberOfRecombinations),
                                              "Number of Extractions = " + str(numberOfExtractions)])
        #print("DURING THIS RUN:")
        #print("Slowest Event Considered =", slowestEvent)
        #print("Fastest Event Considered =", fastestEvent)
        plotEventTimeDistribution(eventLog, outputFiguresDir, fastestEvent, slowestEvent)
        if parameterDict['recordCarrierHistory'] is True:
            plotCarrierTrajectories(allCarriers, parameterDict, deviceArray, outputFiguresDir)
        exit()
    ### THE FOLLOWING CODE WAS USED TO GENERATE THE ANALYSIS DATA USED IN THE MATTYSUMMARIES EXCITON DYNAMICS STUFF
    #print("Plotting Exciton Characteristics")
    #plt.figure()
    #plt.hist(dissExcitonDisp, bins=15)
    #plt.title('Displacement until Dissociation')
    #plt.xlabel('Displacement, nm')
    #plt.ylabel('Frequency, Arb. U.')
    #plt.savefig(outputFiguresDir + 'dissExcitonDisp.pdf')
    #plt.clf()

    #plt.hist(np.array(dissExcitonTime) * 1E9, bins=15
    #plt.title('Time until Dissociation')
    #plt.xlabel('Time, ns')
    #plt.ylabel('Frequency, Arb. U.')
    #plt.savefig(outputFiguresDir + 'dissExcitonTime.pdf')
    #plt.clf()

    #plt.hist(recExcitonDisp, bins=15)
    #plt.title('Displacement until Recombination')
    #plt.xlabel('Displacement, nm')
    #plt.ylabel('Frequency, Arb. U.')
    #plt.savefig(outputFiguresDir + 'recExcitonDisp.pdf')
    #plt.clf()

    #plt.hist(np.array(recExcitonTime) * 1E10, bins=15)
    #plt.title('Time until Recombination')
    #plt.xlabel('Time, ns')
    #plt.ylabel('Frequency, Arb. U.')
    #plt.savefig(outputFiguresDir + 'recExcitonTime.pdf')
    #plt.clf()

    #plt.hist(numberOfHops, bins=15)
    #plt.title('numberOfHops')
    #plt.savefig(outputFiguresDir + 'numberOfExcitonHops.pdf')
    #plt.clf()

    #print("XDE =", numberOfDissociations/parameterDict['minimumNumberOfPhotoinjections'])
    #print("Mean/SD DissDisp =", np.mean(dissExcitonDisp), "+-", np.std(dissExcitonDisp)/np.sqrt(len(dissExcitonDisp)))
    #print("Mean/SD DissTime =", np.mean(dissExcitonTime), "+-", np.std(dissExcitonTime)/np.sqrt(len(dissExcitonTime)))
    #print("Mean/SD RecDisp =", np.mean(recExcitonDisp), "+-", np.std(recExcitonDisp)/np.sqrt(len(recExcitonDisp)))
    #print("Mean/SD RecTime =", np.mean(recExcitonTime), "+-", np.std(recExcitonTime)/np.sqrt(len(recExcitonTime)))

    #plotConnections(excitonPath, chromophoreData.returnChromophoreList([1,8,6]), morphologyData.returnAAMorphology([1,8,6]), outputFiguresDir)

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
    #plt.savefig(outputFiguresDir + 'excitonMSD.pdf')

    #print("Plotting the extraction data as a function of simulation time")
    #plt.figure()
    #plt.plot(convGlobalTime, convExtractions)
    #plt.savefig(outputFiguresDir + 'convergence.pdf')
    #return
    time = T.time() - t0
    print("SIMULATION COMPLETED")
    print("Plotting output graphs...")
    helperFunctions.writeToFile(logFile, ["---=== RESULTS FROM CPU RANK " + sys.argv[2] + " ===---",
                                          "Simulation voltage = " + str(voltageVal),
                                          "Run completed after " + str(KMCIterations) +
                                          " iterations (globalTime = " + str(globalTime) + ") after " +
                                          str(time) + " seconds",
                                          "Number of Photoinjections = " + str(numberOfPhotoinjections),
                                          "Number of Cathode Injections = " + str(numberOfCathodeInjections),
                                          "Number of Anode Injections = " + str(numberOfAnodeInjections),
                                          "Number of Dissociations = " + str(numberOfDissociations),
                                          "Number of Recombinations = " + str(numberOfRecombinations),
                                          "Number of Extractions = " + str(numberOfExtractions)])
    plotEventTimeDistribution(eventLog, outputFiguresDir, fastestEvent, slowestEvent)
    if parameterDict['recordCarrierHistory'] is True:
        plotCarrierTrajectories(allCarriers, parameterDict, deviceArray, outputFiguresDir)


def slurmTimeInS(slurmTime):
    # Expects slurmTime in HH:MM:SS
    splitTime = slurmTime.split(':')
    timeInS = int(splitTime[0]) * 60**2 + int(splitTime[1]) * 60 + int(splitTime[2])
    return timeInS


if __name__ == "__main__":
    global logFile

    KMCDirectory = sys.argv[1]
    CPURank = int(sys.argv[2])
    R.seed(int(sys.argv[3]))
    overwrite = False
    try:
        overwrite = bool(sys.argv[4])
    except:
        pass

    # Get the time limit
    timeLimit = os.getenv('SLURM_TIMELIMIT', None)
    if timeLimit is not None:
        timeLimit = slurmTimeInS(timeLimit)

    jobsFileName = KMCDirectory + '/KMCData_%02d.pickle' % (CPURank)
    deviceDataFileName = KMCDirectory.replace('/KMC', '/code/deviceData.pickle')

    with open(deviceDataFileName, 'rb') as pickleFile:
        [deviceArray, chromophoreData, morphologyData, parameterDict] = pickle.load(pickleFile)
    with open(jobsFileName, 'rb') as pickleFile:
        jobsToRun = pickle.load(pickleFile)
    if parameterDict['outputLogToSTDOUT'] is True:
        print("Redirecting log to standard out.")
        logFile = 'stdout'
    else:
        logFile = KMCDirectory + '/KMClog_' + str(CPURank) + '.log'
        # Reset the log file
        with open(logFile, 'wb+') as logFileHandle:
            pass
    helperFunctions.writeToFile(logFile, ['Found ' + str(len(jobsToRun)) + ' jobs to run:', repr(jobsToRun)])
    #helperFunctions.writeToFile(logFile, ['Found ' + str(len(jobsToRun)) + ' jobs to run.'])
    # Set the affinities for this current process to make sure it's maximising available CPU usage
    currentPID = os.getpid()
    try:
        affinityJob = sp.Popen(['taskset', '-pc', str(CPURank), str(currentPID)], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
        # helperFunctions.writeToFile(logFile, affinityJob[0].split('\n')) #stdOut for affinity set
        # helperFunctions.writeToFile(logFile, affinityJob[1].split('\n')) #stdErr for affinity set
    except OSError:
        helperFunctions.writeToFile(logFile, ["Taskset command not found, skipping setting of processor affinity..."])
        #helperFunctions.writeToFile(logFile, ["Taskset command not found, skipping setting of processor affinity..."])

    # Begin the simulation
    for voltageVal in jobsToRun:
        execute(deviceArray, chromophoreData, morphologyData, parameterDict, voltageVal, timeLimit)
