import numpy as np
import helperFunctions
import copy
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time as T
import random as R
try:
    import mpl_toolkits.mplot3d.axes3d as p3
except ImportError:
    print "Could not import 3D plotting engine, calling the plotMolecule3D function will result in an error!"
    pass



class obtain:
    def __init__(self, morphologyData, moleculeProperties, boxSize, outputDir, moleculeAAIDs):
        self.morphologyData = morphologyData
        self.moleculeProperties = moleculeProperties # Thio COMs are wrapped
        self.boxSize = boxSize
        self.simDims = [[-boxSize[0]/2.0, boxSize[0]/2.0], [-boxSize[1]/2.0, boxSize[1]/2.0], [-boxSize[2]/2.0, boxSize[2]/2.0]]
        print "BoxSize =", self.boxSize, "SimDims =", self.simDims
        self.maximumHoppingDistance = 10.0 # REMEMBER, this is in angstroems not distance units
        self.getChromoPosns()
        self.updateNeighbourList()
        chromophorePairIDs, chromophorePairs = self.obtainHoppingPairs()

        # chromophorePairIDs = [[12, 373], [12, 374], [13, 375], [14, 376], [15, 377], [16, 373], [17, 378], [18, 379], [19, 380]]
        # for hoppingPair in chromophorePairs:
        #     print "PairType =", hoppingPair[0]['periodic'], hoppingPair[1]['periodic']
        #     print "Real ChromoIDs", hoppingPair[0]['realChromoID'], hoppingPair[1]['realChromoID']
        #     print "First lot of atoms ="
        #     print ' '.join(str(ID) for ID in hoppingPair[0]['atomID'])
        #     print "Second lot of atoms ="
        #     print ' '.join(str(ID) for ID in hoppingPair[1]['atomID'])
        #     raw_input('PAUSE')
        # exit()
        helperFunctions.checkORCAFileStructure(outputDir)
        for chromophorePair in chromophorePairs:
            helperFunctions.writeORCAInp(chromophorePair, outputDir)
        print "\n"


    def getChromoPosns(self):
        self.chromophores = {}
        globalChromoNo = -1 # This number includes any periodic chromophores that are not in the original simulation volume
        self.numberOfPeriodicChromophores = 0
        for molNo, molecule in enumerate(self.moleculeProperties):
            for chromoNo, chromophore in enumerate(molecule['morphologyChromophores']):
                # print "Molecule =", molNo, "ChromophoreNumber in mol =", chromoID, "Actual Chromo Number =", globalChromoNo
                globalChromoNo += 1
                chromoDict = {'molID': molNo, 'chromoID': globalChromoNo, 'periodic': False, 'realChromoID': globalChromoNo, 'unwrapped_position': [], 'position': [], 'image': [0, 0, 0], 'type': [], 'mass': [], 'COMPosn': 0, 'thioCOMs': [], 'neighbours': [], 'atomID':[]}
                # For the positions, I need all of the chromophores to be intact, but in the original image (kinda like VMD displays)
                # To do this, calculate the COM of the chromophore based on the unwrapped coordinates, and then wrap everything back into the box
                for atomID in chromophore:
                    chromoDict['atomID'].append(atomID)
                    chromoDict['unwrapped_position'].append(self.morphologyData['unwrapped_position'][atomID])
                    chromoDict['type'].append(self.morphologyData['type'][atomID])
                    chromoDict['mass'].append(self.morphologyData['mass'][atomID])
                # To calculate the COMPosn and avoid cross-boundary issues:
                # Calc it from the unwrapped positions, then wrap it back into the box
                unwrappedThioCOMs = molecule['unwrappedThioCOMs'][chromoNo]
                unwrappedCOMPosn = helperFunctions.calcCOM(chromoDict['unwrapped_position'], chromoDict['mass'])
                imageMultiplier = [0, 0, 0] # To bring it back into the box
                for axis in [0, 1, 2]:
                    while unwrappedCOMPosn[axis] < self.simDims[axis][0]:
                        imageMultiplier[axis] += 1
                        unwrappedCOMPosn[axis] += self.boxSize[axis]
                    while unwrappedCOMPosn[axis] > self.simDims[axis][1]:
                        imageMultiplier[axis] -= 1
                        unwrappedCOMPosn[axis] -= self.boxSize[axis]
                for position in chromoDict['unwrapped_position']:
                    chromoDict['position'].append([position[0] + (imageMultiplier[0]*self.boxSize[0]), position[1] + (imageMultiplier[1]*self.boxSize[1]), position[2] + (imageMultiplier[2]*self.boxSize[2])])
                for thioCOM in unwrappedThioCOMs:
                    chromoDict['thioCOMs'].append([thioCOM[0] + (imageMultiplier[0]*self.boxSize[0]), thioCOM[1] + (imageMultiplier[1]*self.boxSize[1]), thioCOM[2] + (imageMultiplier[2]*self.boxSize[2])])
                chromoDict['COMPosn'] = unwrappedCOMPosn
                self.chromophores[globalChromoNo] = chromoDict
                # Add in a periodic segment if the COM position is within self.maximumHoppingDistance of a boundary
                periodLims = [[self.simDims[0][0]-self.maximumHoppingDistance, self.simDims[0][1]+self.maximumHoppingDistance], [self.simDims[1][0]-self.maximumHoppingDistance, self.simDims[1][1]+self.maximumHoppingDistance], [self.simDims[2][0]-self.maximumHoppingDistance, self.simDims[2][1]+self.maximumHoppingDistance]]
                realChromoID = globalChromoNo
                for ximage in range(-1,2,1):
                    for yimage in range(-1,2,1):
                        for zimage in range(-1,2,1):
                            if (ximage == 0) and (yimage == 0) and (zimage == 0):
                                continue
                            periodicCOMX = chromoDict['COMPosn'][0]+(ximage*self.boxSize[0])
                            periodicCOMY = chromoDict['COMPosn'][1]+(yimage*self.boxSize[1])
                            periodicCOMZ = chromoDict['COMPosn'][2]+(zimage*self.boxSize[2])
                            if ((periodicCOMX >= periodLims[0][0]) and (periodicCOMX <= periodLims[0][1])) and ((periodicCOMY >= periodLims[1][0]) and (periodicCOMY <= periodLims[1][1])) and ((periodicCOMZ >= periodLims[2][0]) and (periodicCOMZ <= periodLims[2][1])):
                                # print "PeriodLims =", periodLims
                                # print "Original COM =", chromoDict['COMPosn']
                                # print "New COM at image", [ximage, yimage, zimage], '=', [periodicCOMX, periodicCOMY, periodicCOMZ]
                                # This chromophore is within the self.maximumHoppingDistance of a boundary so add its periodic partners to the dictionary!
                                globalChromoNo += 1
                                self.numberOfPeriodicChromophores += 1
                                periodicChromoDict = copy.deepcopy(chromoDict)
                                periodicChromoDict['realChromoID'] = realChromoID
                                periodicChromoDict['chromoID'] = globalChromoNo
                                periodicChromoDict['periodic'] = True
                                periodicChromoDict['image'] = [ximage, yimage, zimage]
                                periodicChromoDict['COMPosn'] = np.array([periodicCOMX, periodicCOMY, periodicCOMZ])
                                previousPos = chromoDict['position']
                                for atomID, position in enumerate(periodicChromoDict['position']):
                                    periodicChromoDict['position'][atomID] = [position[0] + (ximage*self.boxSize[0]), position[1] + (yimage*self.boxSize[1]), position[2] + (zimage*self.boxSize[2])]
                                finalPos = periodicChromoDict['position']
                                previousThioID = chromoDict['thioCOMs']
                                for thioID, position in enumerate(periodicChromoDict['thioCOMs']):
                                    periodicChromoDict['thioCOMs'][thioID] = [position[0] + (ximage*self.boxSize[0]), position[1] + (yimage*self.boxSize[1]), position[2] + (zimage*self.boxSize[2])]
                                # finalThioID = periodicChromoDict['thioCOMs']
                                # print [ximage, yimage, zimage]
                                # print "ATOM POSITIONS"
                                # for i in range(len(previousPos)):
                                #     print previousPos[i], "\t\t", finalPos[i]
                                # print "THIO POSITIONS"
                                # for i in range(len(previousThioID)):
                                #     print previousThioID[i], "\t\t", finalThioID[i]
                                # raw_input('HALT')
                                self.chromophores[globalChromoNo] = periodicChromoDict


    def updateNeighbourList(self):
        # Mike's clustering experimentation might help a lot with this, so for now I'm going to just O(N^{2}) brute force it.
        # This now checks the centres of mass of each thiophene group with each other thiophene group to see if any are within
        # self.maximumHoppingDistance of each other.
        print "Determining neighbours for", len(self.chromophores), "chromophores..."
        t0 = T.time()
        for chromophore1ID in self.chromophores.keys():
            for chromophore2ID in self.chromophores.keys():
                if chromophore1ID == chromophore2ID:
                    continue
                neighbourAdded = False
                for thio1 in self.chromophores[chromophore1ID]['thioCOMs']:
                    if neighbourAdded == True:
                       break
                    for thio2 in self.chromophores[chromophore2ID]['thioCOMs']:
                        if neighbourAdded == True:
                           break
                        if helperFunctions.calculateSeparation(thio1, thio2) <= self.maximumHoppingDistance:
                            #print self.chromophores[chromophore1ID]['chromoID'], "and", self.chromophores[chromophore2ID]['chromoID'], 'match with', thio1, '->', thio2, '<', self.maximumHoppingDistance
                            if self.chromophores[chromophore2ID]['chromoID'] not in self.chromophores[chromophore1ID]['neighbours']:
                                self.chromophores[chromophore1ID]['neighbours'].append(self.chromophores[chromophore2ID]['chromoID'])
                            neighbourAdded = True
        t1 = T.time()
        print "Neighbour calculations took %.2f seconds." % (t1-t0)


    def obtainHoppingPairs(self):
        chromophorePairs = []
        chromophorePairIDs = []
        print "Obtaining hopping pairs based on neighbouring chromophores..."
        t0 = T.time()
        for chromoID in self.chromophores:
            for neighbour in self.chromophores[chromoID]['neighbours']:
                if self.chromophores[neighbour]['realChromoID'] == self.chromophores[chromoID]['realChromoID']:
                    print "Chromophore 1 =", self.chromophores[chromoID]
                    print "Chromophore 2 =", self.chromophores[neighbour]
                    print "ChromoIDs =", self.chromophores[chromoID]['chromoID'], self.chromophores[neighbour]['chromoID']
                    print "Neighbours =", self.chromophores[chromoID]['neighbours']
                    raise SystemError("Hopping from one chromo to its own periodic image")
                else:
                    # Don't need to include the backwards hop
                    if chromoID < neighbour:
                        if [chromoID, neighbour] not in chromophorePairIDs:
                            chromophorePairIDs.append([chromoID, neighbour])
                    else:
                        if [neighbour, chromoID] not in chromophorePairIDs:
                            chromophorePairIDs.append([neighbour, chromoID])
        chromophorePairIDs.sort()
        # Many of these hop pairs are duplicates, but they all give the same tranfer integral, so we only need to have one of them in the hop pair list.
        realIDs = {}
        for hopPair in chromophorePairIDs:
            chromo1ID = hopPair[0]
            chromo2ID = hopPair[1]
            realChromo1ID = self.chromophores[chromo1ID]['realChromoID']
            realChromo2ID = self.chromophores[chromo2ID]['realChromoID']
            dictKey = str([realChromo1ID, realChromo2ID])
            if dictKey not in realIDs:
                realIDs[dictKey] = [[chromo1ID, chromo2ID]]
            else:
                realIDs[dictKey].append([chromo1ID, chromo2ID])
        for key in realIDs.keys():
            if len(realIDs[key]) > 1:
                for item in realIDs[key][1:]:
                    chromophorePairIDs.remove(item)
        # This leaves the first transfer integral for the hop but ignores the others
        for chromophorePair in chromophorePairIDs:
            chromophorePairs.append([self.chromophores[chromophorePair[0]], self.chromophores[chromophorePair[1]]])
        t1 = T.time()
        print "Pairing treatment took %.2f seconds." % (t1-t0)
        print "There are", len(chromophorePairIDs), "pairs of chromophores"
        return chromophorePairIDs, chromophorePairs


def plotMolecule3D(chromophores, simDims, chromoNumbers=None):
    chromoLocs = []
    thioLocs = []
    if chromoNumbers != None:
        thioLocs1 = []
        thioLocs2 = []
        atomLocs1 = []
        atomLocs2 = []
        chromo1 = chromoNumbers[0]
        chromo2 = chromoNumbers[1]
        for thioPosn in chromophores[chromo1]['thioCOMs']:
            thioLocs1.append(thioPosn)
        for thioPosn in chromophores[chromo2]['thioCOMs']:
            thioLocs2.append(thioPosn)
        # for atomPosn in chromophores[chromo1]['position']:
        #     atomLocs1.append(atomPosn)
        # for atomPosn in chromophores[chromo2]['position']:
        #     atomLocs2.append(atomPosn)
        print "\n"
        print "THIS GRAPH IS", chromo1, chromo2
    # for chromophore in chromophores.values():
    #     chromoLocs.append(chromophore['COMPosn'])
    #     thioLocs.append(chromophore['thioCOMs'])
    # chromoLocs = np.array(chromoLocs)
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    for thiophene in thioLocs1:
        ax.scatter(thiophene[0], thiophene[1], thiophene[2], s = 50, c = 'r')
    for thiophene in thioLocs2:
        ax.scatter(thiophene[0], thiophene[1], thiophene[2], s = 50, c = 'b')
    # for atom in atomLocs1:
    #     ax.scatter(atom[0], atom[1], atom[2], s = 20, c = 'g')
    # for atom in atomLocs2:
    #     ax.scatter(atom[0], atom[1], atom[2], s = 20, c = 'k')
    #ax.scatter(allAtoms[:,0], allAtoms[:,1], allAtoms[:,2], s = 20, c = 'g')
    #ax.scatter(chromoLocs[:,0], chromoLocs[:,1], chromoLocs[:,2], s = 50, c = 'r')
    # for chromophore in thioLocs:
    #     chromox = []
    #     chromoy = []
    #     chromoz = []
    #     for thiophene in chromophore:
    #         chromox.append(thiophene[0])
    #         chromoy.append(thiophene[1])
    #         chromoz.append(thiophene[2])
    #     rval = R.random()
    #     gval = R.random()
    #     bval = R.random()
    #     ax.scatter(chromox, chromoy, chromoz, s = 50, c = (rval, gval, bval))
    #ax.set_xlim((-8,8))
    ax.set_xlim((simDims[0][0], simDims[0][1]))
    ax.set_ylim((simDims[1][0], simDims[1][1]))
    #ax.set_zlim((-15,15))
    ax.set_zlim((simDims[2][0], simDims[2][1]))
    print "Wobbey show"
    plt.show()
