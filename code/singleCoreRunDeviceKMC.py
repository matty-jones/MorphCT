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
        self.prefactor = parameterDict['hoppingPrefactor']
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
            hopRate = calculateFRETHopRate(self.prefactor, self.lifetimeParameter, self.rF, rij, deltaEij * elementaryCharge, self.T)
            hopTime = determineHopTime(hopRate)
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
                print("Out of bounds Exciton hop discarded")
            else:
                self.currentDevicePosn = list(np.array(self.currentDevicePosn) + np.array(self.destinationImage))
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
        self.prefactor = parameterDict['hoppingPrefactor']
        self.recombining = False
        self.recombiningWith = None
        self.T = parameterDict['systemTemperature']
        # self.carrierType: Set carrier type to be == 0 if Electron and == 1 if Hole. This allows us to do quick arithmetic to get the signs correct in the potential calculations without having to burn through a ton of conditionals.
        if self.currentChromophore.species == 'Donor':
            self.carrierType = HOLE
            self.lambdaij = parameterDict['reorganisationEnergyDonor']
        elif self.currentChromophore.species == 'Acceptor':
            self.carrierType = ELECTRON
            self.lambdaij = parameterDict['reorganisationEnergyAcceptor']
        self.relativePermittivity = parameterDict['relativePermittivity']

    def calculateBehaviour(self):
        try:
            [self.destinationChromophore, self.hopTime, self.destinationImage] = self.calculateHop()
        except ValueError:
            [self.destinationChromophore, self.hopTime, self.destinationImage] = [None, None, None]

    def calculateInjectHop(self, destinationChromophore, deltaE, transferIntegral):
        self.destinationChromophore = copy.deepcopy(destinationChromophore)
        hoppingRate = calculateMarcusHopRate(self.prefactor, self.lambdaij * elementaryCharge, transferIntegral * elementaryCharge, deltaE, self.T)
        self.hopTime = determineHopTime(hoppingRate)
        self.destinationImage = [0, 0, 0]

    def calculateHop(self):
        # Determine the hop times to all possible neighbours
        hopTimes = []
        # Obtain the reorganisation energy in J (from eV in the parameter file)
        for neighbourIndex, transferIntegral in enumerate(self.currentChromophore.neighboursTI):
            # Ignore any hops with a NoneType transfer integral (usually due to an ORCA error)
            if transferIntegral is None:
                continue
            destinationChromophore = globalChromophoreData.returnChromophoreList(self.currentDevicePosn)[self.currentChromophore.neighbours[neighbourIndex][0]]
            # Need to determine the relative image between the two (recorded in obtainChromophores) to check if the carrier is hopping out of this morphology cell and into an adjacent one
            neighbourRelativeImage = self.currentChromophore.neighbours[neighbourIndex][1]
            deltaEij = self.calculatePotential(destinationChromophore, neighbourRelativeImage, self.currentChromophore.neighboursDeltaE[neighbourIndex])
            # All of the energies (EXCEPT EIJ WHICH IS ALREADY IN J) are in eV currently, so convert them to J
            hopRate = calculateMarcusHopRate(self.prefactor, self.lambdaij * elementaryCharge, transferIntegral * elementaryCharge, deltaEij, self.T)
            hopTime = determineHopTime(hopRate)
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
            # We're hopping over a boundary. If the closest chromophore to the destination position in the adjacent device cell is the same type as the one we're currently on, permit the hop with the already-calculated hopTime. Otherwise, ignore this hop.
            targetDevicePosn = list(np.array(self.currentDevicePosn) + np.array(self.destinationImage))
            newChromophore = globalChromophoreData.returnClosestChromophoreToPosition(targetDevicePosn, destinationPosition)
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
                if self.carrierType == ELECTRON:
                    print('EVENT: Electron Left Device! New number of extractions:', numberOfExtractions)
                else:
                    print('EVENT: Hole Left Device! New number of extractions:', numberOfExtractions)
            elif newChromophore == 'Out of Bounds':
                # Trying to cross a periodic boundary. Das ist verboten, discarding the hop
                pass
            elif newChromophore.species == self.currentChromophore.species:
                self.currentDevicePosn = list(np.array(self.currentDevicePosn) + np.array(self.destinationImage))
                self.currentChromophore = newChromophore
            #else:
            #    print("Crossing a boundary to a different chromophore species. Das ist verboten, discarding the hop")

    def calculatePotential(self, destinationChromophore, neighbourRelativeImage, chromoEij):
        global globalCarrierDict

        # The potential will be calculated in J
        currentAbsolutePosition = np.array(self.currentDevicePosn)*parameterDict['morphologyCellSize'] + (np.array(self.currentChromophore.posn) * 1E-10)
        destinationZ = ((np.array(self.currentDevicePosn) + np.array(neighbourRelativeImage))*parameterDict['morphologyCellSize'] + (np.array(destinationChromophore.posn) * 1E-10))[2]
        zSep = destinationZ - currentAbsolutePosition[2]
        charge = elementaryCharge * ((2 * self.carrierType) - 1)
        totalPotential = (zSep * currentFieldValue * charge) + (chromoEij * charge)
        # Coulombic component!
        coulombicPotential = 0.0
        coulombConstant = 1.0 / (4 * np.pi * epsilonNought * self.relativePermittivity)
        for carrier in globalCarrierDict.values():
            # Each carrier feels coulombic effects from every other carrier in the system, along with the image charges induced in the top and bottom contacts, AND the image charges induced by the image charges.
            # The image charges are needed in order to keep the contacts as perfect metals (i.e. no internal electric field)
            # Check papers from Chris Groves and Ben Lyons for more details.
            # TODO I can't work out how image charges work (and whether their effect is negligible anyway) so for now I'm just ignoring them
            carrierPosn = np.array(carrier.currentDevicePosn)*parameterDict['morphologyCellSize'] + (np.array(carrier.currentChromophore.posn) * 1E-10)
            if carrier.ID != self.ID:
                # Only consider the current carrier if it is not us!
                separation = np.linalg.norm(currentAbsolutePosition - carrierPosn)
                # TODO If separation < 1E-9, create a recombination event here in the carriers
                coulombicPotential += coulombConstant * ((elementaryCharge * ((2 * self.carrierType) - 1)) * (elementaryCharge * ((2 * carrier.carrierType) - 1)) / separation)
                # I'm also going to use this opportunity (as we're iterating over all carriers in the system) to see if we're close enough to any to recombine.
                if self.recombining is False:
                    # Only do this if we're not currently recombining with something
                    if (separation <= parameterDict['coulombCaptureRadius']) and (self.carrierType != carrier.carrierType):
                        # If carriers are within 1nm of each other, then assume that they are within the Coulomb capture radius and are about to recombine
                         self.recombining = True
                         self.recombiningWith = carrier.ID
                         # We don't necessarily need to update the carrier.recombining, since all of the behaviour will be determined in the main program by checking hoppingCarrier.recombining, but this is more technically correct.
                         carrier.recombining = True
                         carrier.recombiningWith = self.ID
            ## Now consider the image charges and the images of the images.
            ## First do the top image charge
            #deviceZSize = globalChromophoreData.deviceArray.shape[2] * parameterDict['morphologyCellSize']
            #topImagePosn = copy.deepcopy(carrierPosn)
            #topImagePosn[2] = deviceZSize + (deviceZSize - carrierPosn[2])
            #separation = np.linalg.norm(currentAbsolutePosition - topImagePosn)
            #coulombicPotential -= ((2 * self.carrierType) - 1) * (coulombConstant / separation)
            ## Now do the bottom image charge
            #bottomImagePosn = copy.deepcopy(carrierPosn)
            #bottomImagePosn[2] = - carrierPosn[2]
            #separation = np.linalg.norm(currentAbsolutePosition - bottomImagePosn)
            #coulombicPotential -= ((2 * self.carrierType) - 1) * (coulombConstant / separation)
            ## Now do the top image of the bottom image charge
            #topImageOfBottomImagePosn = copy.deepcopy(carrierPosn)
            #topImageOfBottomImagePosn[2] = (2 * deviceZSize) + carrierPosn[2]
            #separation = np.linalg.norm(currentAbsolutePosition - topImageOfBottomImagePosn)
            #coulombicPotential -= ((2 * self.carrierType) -1 ) * (coulombConstant / separation)
            ## And finally the bottom image of the top image charge
            #bottomImageOfTopImagePosn = copy.deepcopy(carrierPosn)
            #bottomImageOfTopImagePosn[2] = - (deviceZSize + (deviceZSize - carrierPosn[2]))
            #separation = np.linalg.norm(currentAbsolutePosition - bottomImageOfTopImagePosn)
            #coulombicPotential -= ((2 * self.carrierType) - 1) * (coulombConstant / separation)
        #print(currentAbsolutePosition)
        #print(currentAbsolutePosition[2], destinationZ, zSep)
        #print(coulombicPotential, zSep * currentFieldValue * charge, chromoEij * charge)
        totalPotential += coulombicPotential
        return totalPotential


class injectSite:
    def __init__(self, devicePosn, chromophore, injectRate, electrode):
        self.devicePosn = devicePosn
        self.chromophore = chromophore
        self.injectRate = injectRate
        self.electrode = electrode
        self.injectTime = self.calculateInjectTime()

    def calculateInjectTime(self):
        return determineHopTime(self.injectRate)


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


def calculateDarkCurrentInjectRates(deviceArray, parameterDict):
    # NOTE In the test bilayer device, there is no way to inject holes from the cathode or electrons from the anode
    # (which is the `easy' injection pathway). As such, there will be no dark current injection in the system
    injectQueue = []
    # Get the device shape so that we know which cells to inject into
    deviceShape = deviceArray.shape
    # Calculate the important energy levels
    bandgap = parameterDict['acceptorLUMO'] - parameterDict['donorHOMO']
    electronInjectBarrier = parameterDict['acceptorLUMO'] - parameterDict['cathodeWorkFunction']
    holeInjectBarrier = parameterDict['anodeWorkFunction'] - parameterDict['donorHOMO']
    # Consider bottom electrode (cathode) first:
    zVal = 0
    for xVal in range(deviceShape[0]):
        for yVal in range(deviceShape[1]):
            injectChromophores = []
            AAMorphology = globalMorphologyData.returnAAMorphology([xVal, yVal, zVal])
            morphologyChromophores = globalChromophoreData.returnChromophoreList([xVal, yVal, zVal])
            estimatedTransferIntegral = estimateTransferIntegral(morphologyChromophores) * elementaryCharge
            for chromophore in morphologyChromophores:
                # Find all chromophores that are within the bottom nanometer (10 Ang) of the device cell
                if chromophore.posn[2] <= -(AAMorphology['lz'] / 2.0) + 10:
                    injectChromophores.append(chromophore)
            for chromophore in injectChromophores:
                if chromophore.species == 'Donor':
                    # Injecting a hole from the cathode (easy)
                    #print("EASY INJECT", holeInjectBarrier, "+", chromophore.HOMO - parameterDict['donorHOMO'])
                    #print("Chromo HOMO =", chromophore.HOMO)
                    #print("Target HOMO =", parameterDict['donorHOMO'])
                    deltaE = elementaryCharge * (holeInjectBarrier + (chromophore.HOMO - parameterDict['donorHOMO']))
                    #print("Easy HOLE Inject from Cathode =", deltaE / elementaryCharge)
                    #input("PAUSE...")
                else:
                    # Injecting an electron from the cathode (hard)
                    #print("Difficult INJECT", bandgap, "-", holeInjectBarrier, "-", (chromophore.HOMO - parameterDict['donorHOMO']))
                    #print("Chromo HOMO =", chromophore.HOMO)
                    #print("Target HOMO =", parameterDict['donorHOMO'])
                    # Here, we make the assumption that the energy difference between the chromophore and literature HOMOs
                    # would be the same as the energy difference between the chromophore and literature LUMOs for injection
                    # (this is because we don't record the Donor LUMO or the Acceptor HOMO)
                    deltaE = elementaryCharge * (bandgap - holeInjectBarrier - (chromophore.HOMO - parameterDict['donorHOMO']))
                    #print("Difficult HOLE inject from Cathode =", deltaE / elementaryCharge)
                    #input("PAUSE...")
                injectRate = calculateMarcusHopRate(parameterDict['hoppingPrefactor'], parameterDict['reorganisationEnergyDonor'] * elementaryCharge, estimatedTransferIntegral, deltaE, parameterDict['systemTemperature'])
                # To prevent a queue filled with useless difficult injection events that will never happen, only add rates > 1
                if injectRate >= 1e0:
                    injection = injectSite([xVal, yVal, zVal], chromophore, injectRate, 'Cathode')
                    heapq.heappush(injectQueue, (injection.injectTime, injection))
    # Now consider top electrode (anode):
    zVal = deviceShape[2]-1
    for xVal in range(deviceShape[0]):
        for yVal in range(deviceShape[1]):
            injectChromophores = []
            AAMorphology = globalMorphologyData.returnAAMorphology([xVal, yVal, zVal])
            morphologyChromophores = globalChromophoreData.returnChromophoreList([xVal, yVal, zVal])
            estimatedTransferIntegral = estimateTransferIntegral(morphologyChromophores) * elementaryCharge
            for chromophore in morphologyChromophores:
                # Find all chromophores that are within the top nanometer (10 Ang) of the device cell
                if chromophore.posn[2] >= (AAMorphology['lz'] / 2.0) - 10:
                    injectChromophores.append(chromophore)
            for chromophore in injectChromophores:
                if chromophore.species == 'Donor':
                    # Injecting a hole from the anode (hard)
                    deltaE = elementaryCharge * (-electronInjectBarrier + bandgap + (chromophore.LUMO - parameterDict['acceptorLUMO']))
                else:
                    # Injecting an electron from the anode (easy)
                    deltaE = elementaryCharge * (electronInjectBarrier - (chromophore.LUMO - parameterDict['acceptorLUMO']))
                injectRate = calculateMarcusHopRate(parameterDict['hoppingPrefactor'], parameterDict['reorganisationEnergyAcceptor'] * elementaryCharge, estimatedTransferIntegral, deltaE, parameterDict['systemTemperature'])
                # To prevent a queue filled with useless difficult injection events that will never happen, only add rates > 1
                if injectRate >= 1e0:
                    injection = injectSite([xVal, yVal, zVal], chromophore, injectRate, 'Anode')
                    heapq.heappush(injectQueue, (injection.injectTime, injection))
    return injectQueue


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
    #    separation = np.linalg.norm(currentPosn - chromophoreList[chromoID].posn)
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


def calculateMarcusHopRate(prefactor, lambdaij, Tij, deltaEij, T):
    ## Semiclassical Marcus Hopping Rate Equation
    #print("Lambdaij =", lambdaij)
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


def determineHopTime(rate):
    # Use the KMC algorithm to determine the wait time to this hop
    if rate != 0:
        while True:
            x = R.random()
            # Ensure that we do not get exactly 0.0 or 1.0, which would break our logarithm
            if (x != 0.0) and (x != 1.0):
                break
        tau = - np.log(x) / rate
    else:
        # If rate == 0, then make the hopping time extremely long
        tau = 1E99
    return tau


def execute(deviceArray, chromophoreData, morphologyData, parameterDict, voltageVal):
    # ---=== PROGRAMMER'S NOTE ===---
    # This `High Resolution' version of the code will permit chromophore-based hopping through the device, rather than just approximating the distribution. We have the proper resolution there, let's just use it!
    # ---=========================---
    global globalChromophoreData
    global globalMorphologyData
    global globalTime
    global currentFieldValue

    globalChromophoreData = chromophoreData
    globalMorphologyData = morphologyData
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

    # DEBUG
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
    numberOfDarkinjections = 0
    excitonIndex = 0
    carrierIndex = 0
    outputCurrentConverged = False
    recombiningCarrierIDs = []

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
    KMCIterations = 0


    # As the morphology is not changing, we can calculate the dark inject rates at the beginning and then not have to do it again
    darkInjectQueue = calculateDarkCurrentInjectRates(deviceArray, parameterDict)


    print("\n\n ---=== MAIN KMC LOOP START ===---\n\n")
    # Main KMC loop
    while True:
        KMCIterations += 1
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
            # The darkInjectQueue is already ordered based on the most likely injection to occur next.
            # Queue up this dark injection in the main KMC queue. However, because there's not anything stopping us from
            # reinjecting to this site later on, recalculate a new injection time for this site and queue it back up
            # in the darkInjectQueue
            # In some cases (like the test morphology), carriers can't be injected with a rate > 1, and so the darkInjectQueue
            # is empty. This catches that and effectively deactivates dark current if it can't happen.
            if len(darkInjectQueue) > 0:
                nextDarkInject = darkInjectQueue.pop(0)
                # Decrement the dark inject queue
                decrementTime(darkInjectQueue, nextDarkInject.injectTime)
                # Now requeue this inject site
                newDarkInject = copy.deepcopy(nextDarkInject).calculateInjectTime()
                heap.heappush(darkInjectQueue, (newDarkInject.injectTime, newDarkInject))
                # And finally, push this injection to the main KMC queue
                heapq.heappush(eventQueue, (nextDarkInject.injectTime, 'dark', nextDarkInject))

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
            print("EVENT: Photoinjection #" + str(numberOfPhotoinjections), "into", randomDevicePosition, "(which has type", repr(deviceArray[tuple(randomDevicePosition)]) + ")", "after", KMCIterations, "iterations")
            injectedExciton = exciton(excitonIndex, globalTime, randomDevicePosition, parameterDict)
            if (injectedExciton.canDissociate is True) or (injectedExciton.hopTime is None):
                # Injected onto either a dissociation site or a trap site
                if injectedExciton.canDissociate is True:
                    print("EVENT: Exciton Dissociating immediately")
                    numberOfDissociations += 1
                    numberOfHops.append(injectedExciton.numberOfHops)
                    # Create the carrier instances, but don't yet calculate their behaviour (need to add them to the carrier list before we can calculate the energetics)
                    # Also add the carriers to the carrier dictionary for when we need to calc deltaE in the device
                    injectedElectron = carrier(carrierIndex, globalTime, injectedExciton.currentDevicePosn, injectedExciton.electronChromophore, 'Exciton', parameterDict)
                    globalCarrierDict[carrierIndex] = injectedElectron
                    carrierIndex += 1
                    injectedHole = carrier(carrierIndex, globalTime, injectedExciton.currentDevicePosn, injectedExciton.holeChromophore, 'Exciton', parameterDict)
                    globalCarrierDict[carrierIndex] = injectedHole
                    carrierIndex += 1
                    # Now determine the behaviour of both carriers, and add their next hops to the KMC queue
                    injectedElectron.calculateBehaviour()
                    if injectedElectron.hopTime is not None:
                        heapq.heappush(eventQueue, (injectedElectron.hopTime, 'carrierHop', injectedElectron))
                    injectedHole.calculateBehaviour()
                    if injectedHole.hopTime is not None:
                        heapq.heappush(eventQueue, (injectedHole.hopTime, 'carrierHop', injectedHole))
                else:
                    # Injected onto a site with no connections, so this exciton will eventually die
                    print("EVENT: Exciton Recombining immediately")
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
            heapq.heappush(eventQueue, (photoinjectionTime, 'photo', None))
            numberOfPhotoinjections += 1

        elif nextEvent[1] == 'dark':
            injectSite = nextEvent[2]
            print("EVENT: Dark Current injection #" + str(numberOfDarkinjections), "into", injectSite.devicePosn, "(which has type", repr(deviceArray[tuple(injectSite.devicePosn)]) + ")", "after", KMCIterations, "iterations")
            if injectSite.chromophore.species == 'Donor':
                # Inject a hole
                injectedCarrier = carrier(carrierIndex, globalTime, injectSite.devicePosn, injectSite.chromophore, injectSite.electrode, parameterDict)
                globalCarrierDict[carrierIndex] = injectedCarrier
                carrierIndex += 1
            else:
                # Inject an electron
                injectedCarrier = carrier(carrierIndex, globalTime, injectSite.devicePosn, injectSite.chromophore, injectSite.electrode, parameterDict)
                globalCarrierDict[carrierIndex] = injectedCarrier
                carrierIndex += 1
            # Requeue the inject site
            newDarkInject = copy.deepcopy(injectSite).calculateInjectTime()
            heap.heappush(darkInjectQueue, (newDarkInject.injectTime, newDarkInject))
            # Determine the carrier's next hop and queue it
            injectedCarrier.calculateBehaviour()
            if injectedCarrier.hopTime is not None:
                heapq.heappush(eventQueue, (injectedCarrier.hopTime, 'carrierHop', injectedCarrier))
            numberOfDarkinjections += 1

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
                    print("EVENT: Exciton Dissociating", "after", KMCIterations, "iterations")
                    numberOfDissociations += 1
                    numberOfHops.append(injectedExciton.numberOfHops)
                    # Create the carrier instances, but don't yet calculate their behaviour (need to add them to the carrier list before we can calculate the energetics)
                    # Also add the carriers to the carrier dictionary for when we need to calc deltaE in the device
                    # Plop an electron down on the electronChromophore of the dissociating exciton, and a hole on the holeChromophore
                    injectedElectron = carrier(carrierIndex, globalTime, hoppingExciton.currentDevicePosn, hoppingExciton.electronChromophore, 'Exciton', parameterDict)
                    globalCarrierDict[carrierIndex] = injectedElectron
                    carrierIndex += 1
                    injectedHole = carrier(carrierIndex, globalTime, hoppingExciton.currentDevicePosn, hoppingExciton.holeChromophore, 'Exciton', parameterDict)
                    globalCarrierDict[carrierIndex] = injectedHole
                    carrierIndex += 1
                    # Now determine the behaviour of both carriers, and add their next hops to the KMC queue
                    injectedElectron.calculateBehaviour()
                    if injectedElectron.hopTime is not None:
                        heapq.heappush(eventQueue, (injectedElectron.hopTime, 'carrierHop', injectedElectron))
                    injectedHole.calculateBehaviour()
                    if injectedHole.hopTime is not None:
                        heapq.heappush(eventQueue, (injectedHole.hopTime, 'carrierHop', injectedHole))
                else:
                    print("EVENT: Exciton Recombining", "after", KMCIterations, "iterations")
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
            else:
                # Normal carrier hop, so requeue this carrier
                heapq.heappush(eventQueue, (hoppingCarrier.hopTime, 'carrierHop', hoppingCarrier))
            # Check if we're eligible to recombine
            if (hoppingCarrier.recombining is True) and (hoppingCarrier.ID not in recombiningCarrierIDs) and (hoppingCarrier.recombiningWith not in recombiningCarrierIDs):
                recombiningCarrierIDs.append(hoppingCarrier.ID)
                recombiningCarrierIDs.append(hoppingCarrier.recombiningWith)
                recombinationTime = -np.log(R.random()) / parameterDict['recombinationRate']
                heapq.heappush(eventQueue, (recombinationTime, 'recombine', hoppingCarrier))

        elif nextEvent[1] == 'recombine':
            print(eventQueue)
            print("EVENT: Carrier Recombination Check", "after", KMCIterations, "iterations")
            # A recombination event is about to occur. At this point, we should check if the carrier and its recombination partner are still in range.
            carrier1 = nextEvent[2]
            carrier1Posn = np.array(carrier1.currentDevicePosn) * parameterDict['morphologyCellSize'] + (np.array(carrier1.currentChromophore.posn) * 1E-10)
            carrier2 = globalCarrierDict[carrier1.recombiningWith]
            carrier2Posn = np.array(carrier2.currentDevicePosn) * parameterDict['morphologyCellSize'] + (np.array(carrier2.currentChromophore.posn) * 1E-10)
            separation = np.linalg.norm(carrier2Posn - carrier1Posn)
            recombiningCarrierIDs.remove(carrier1.ID)
            recombiningCarrierIDs.remove(carrier2.ID)
            if separation <= parameterDict['coulombCaptureRadius']:
                print(separation, "<=", parameterDict['coulombCaptureRadius'])
                print("EVENT: Carrier Recombination Succeeded", "after", KMCIterations, "iterations")
                # Carriers are in range, so recombine them
                carrier1.removedTime = globalTime
                carrier2.removedTime = globalTime
                numberOfRecombinations += 1
                globalCarrierDict.pop(carrier1.ID)
                globalCarrierDict.pop(carrier2.ID)
            else:
                print(separation, ">", parameterDict['coulombCaptureRadius'])
                print("EVENT: Carrier Recombination Failed", "after", KMCIterations, "iterations")
                # Carriers are no longer in range, so the recombination fails. Update their recombination flags
                carrier1.recombining = False
                carrier1.recombiningWith = None
                carrier2.recombining = False
                carrier2.recombiningWith = None

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
    for voltageVal in jobsToRun:
        execute(deviceArray, chromophoreData, morphologyData, parameterDict, voltageVal)
