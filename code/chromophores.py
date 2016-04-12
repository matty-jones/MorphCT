import numpy as np
import helperFunctions
import copy
import matplotlib.pyplot as plt
import time as T
try:
    import mpl_toolkits.mplot3d.axes3d as p3
except ImportError:
    print "Could not import 3D plotting engine, calling the plotMolecule3D function will result in an error!"
    pass



class obtain:
    def __init__(self, morphologyData, moleculeProperties, boxSize, outputDir, moleculeAAIDs):
        self.morphologyData = morphologyData
        self.moleculeProperties = moleculeProperties
        self.boxSize = boxSize
        self.simDims = [[-boxSize[0]/2.0, boxSize[0]/2.0], [-boxSize[1]/2.0, boxSize[1]/2.0], [-boxSize[2]/2.0, boxSize[2]/2.0]]
        self.maximumHoppingDistance = 10.0 # REMEMBER, this is in angstroems not distance units
        self.getChromoPosns()
        self.updateNeighbourList()
        chromophorePairs = self.obtainHoppingPairs()
        helperFunctions.checkORCAFileStructure(outputDir)
        print "Centre the dimer pairs by their COM before writing the ORCA file"
        exit()
        for chromophorePair in chromophorePairs:
            helperFunctions.writeORCAInp(chromophorePair, outputDir)
            exit()


    def getChromoPosns(self):
        self.chromophores = {}
        globalChromoNo = -1 # This number includes any periodic chromophores that are not in the original simulation volume
        self.numberOfPeriodicChromophores = 0
        for molNo, molecule in enumerate(self.moleculeProperties):
            for chromoNo, chromophore in enumerate(molecule['morphologyChromophores']):
                # print "Molecule =", molNo, "ChromophoreNumber in mol =", chromoID, "Actual Chromo Number =", globalChromoNo
                globalChromoNo += 1
                chromoDict = {'molID': molNo, 'chromoID': globalChromoNo, 'periodic': False, 'realChromoID': globalChromoNo, 'position': [], 'image': [0, 0, 0], 'type': [], 'mass': [], 'COMPosn': 0, 'thioCOMs': [], 'neighbours': [], 'atomID':[]}
                for atomID in chromophore:
                    chromoDict['atomID'].append(atomID)
                    chromoDict['position'].append(self.morphologyData['position'][atomID])
                    chromoDict['type'].append(self.morphologyData['type'][atomID])
                    chromoDict['mass'].append(self.morphologyData['mass'][atomID])
                chromoDict['COMPosn'] = helperFunctions.calcCOM(chromoDict['position'], chromoDict['mass'])
                chromoDict['thioCOMs'] = molecule['thioCOMs'][chromoNo]
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
                                # This chromophore is within the self.maximumHoppingDistance of a boundary so add its periodic partners to the dictionary!
                                globalChromoNo += 1
                                self.numberOfPeriodicChromophores += 1
                                periodicChromoDict = copy.deepcopy(chromoDict)
                                periodicChromoDict['realChromoID'] = realChromoID
                                periodicChromoDict['chromoID'] = globalChromoNo
                                periodicChromoDict['periodic'] = True
                                periodicChromoDict['image'] = [ximage, yimage, zimage]
                                periodicChromoDict['COMPosn'] = np.array([periodicCOMX, periodicCOMY, periodicCOMZ])
                                for atomID, position in enumerate(periodicChromoDict['position']):
                                    periodicChromoDict['position'][atomID] = [position[0] + (ximage*self.boxSize[0]), position[1] + (yimage*self.boxSize[1]), position[2] + (zimage*self.boxSize[2])]
                                for thioID, position in enumerate(periodicChromoDict['thioCOMs']):
                                    periodicChromoDict['thioCOMs'][thioID] = [position[0] + (ximage*self.boxSize[0]), position[1] + (yimage*self.boxSize[1]), position[2] + (zimage*self.boxSize[2])]
                                self.chromophores[globalChromoNo] = periodicChromoDict


    def updateNeighbourList(self):
        # Mike's clustering experimentation might help a lot with this, so for now I'm going to just O(N^{2}) brute force it.
        # This now checks the centres of mass of each thiophene group with each other thiophene group to see if any are within
        # self.maximumHoppingDistance of each other.
        print "Determining neighbours for", len(self.chromophores), "chromophores..."
        t0 = T.time()
        for chromophore1ID in self.chromophores:
            for chromophore2ID in self.chromophores:
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
                            if self.chromophores[chromophore1ID]['chromoID'] not in self.chromophores[chromophore2ID]['neighbours']:
                                self.chromophores[chromophore2ID]['neighbours'].append(self.chromophores[chromophore1ID]['chromoID'])
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
                elif chromoID < neighbour:
                    if [chromoID, neighbour] not in chromophorePairIDs:
                        chromophorePairIDs.append([chromoID, neighbour])
                else:
                    if [neighbour, chromoID] not in chromophorePairIDs:
                        chromophorePairIDs.append([neighbour, chromoID])
        for chromophorePair in chromophorePairIDs:
            chromophorePairs.append([self.chromophores[chromophorePair[0]], self.chromophores[chromophorePair[1]]])
        t1 = T.time()
        print "Pairing treatment took %.2f seconds." % (t1-t0)
        print "There are", len(chromophorePairIDs), "pairs of chromophores"
        return chromophorePairs







def plotMolecule3D(chromophores, simDims):
    chromoLocs = []
    for chromophore in chromophores.values():
        chromoLocs.append(chromophore['COMPosn'])
    chromoLocs = np.array(chromoLocs)
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    #ax.scatter(allAtoms[:,0], allAtoms[:,1], allAtoms[:,2], s = 20, c = 'g')
    ax.scatter(chromoLocs[:,0], chromoLocs[:,1], chromoLocs[:,2], s = 50, c = 'r')
    #ax.set_xlim((-8,8))
    ax.set_xlim((simDims[0][0], simDims[0][1]))
    ax.set_ylim((simDims[1][0], simDims[1][1]))
    #ax.set_zlim((-15,15))
    ax.set_zlim((simDims[2][0], simDims[2][1]))
    plt.show()
