import numpy as np
import helperFunctions
import copy
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
        self.maximumHoppingDistance = 10.0
        self.getChromoPosns()
        print "CHROMO 0000"
        print self.chromophores[0]
        print "CHROMO 0001"
        print self.chromophores[1]
        exit()
        self.updateNeighbourList()
        chromophorePairs = self.obtainHoppingPairs()
        helperFunctions.checkORCAFileStructure(outputDir)
        for chromophorePair in chromophorePairs:
            helperFunctions.writeORCAInp(chromophorePair, outputDir)
            exit()


    def getChromoPosns(self):
        self.chromophores = {}
        globalChromoNo = -1 # This number includes any periodic chromophores that are not in the original simulation volume
        self.numberOfPeriodicChromophores = 0
        for molNo, molecule in enumerate(self.moleculeProperties):
            for chromoNo, chromophore in enumerate(molecule['morphologyChromophores']):
                globalChromoNo += 1
                chromoDict = {'molID': molNo, 'chromoID': globalChromoNo, 'periodic': False, 'realChromoID': globalChromoNo, 'position': [], 'image': 0, 'type': [], 'mass': [], 'COMPosn': 0, 'neighbours': [], 'atomID':[]}
                for atomID in chromophore:
                    chromoDict['atomID'].append(atomID)
                    chromoDict['position'].append(self.morphologyData['position'][atomID])
                    chromoDict['type'].append(self.morphologyData['type'][atomID])
                    chromoDict['mass'].append(self.morphologyData['mass'][atomID])
                chromoDict['COMPosn'] = helperFunctions.calcCOM(chromoDict['position'], chromoDict['mass'])
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
                                self.chromophores[globalChromoNo] = periodicChromoDict


    def updateNeighbourList(self):
        # Mike's clustering experimentation might help a lot with this, so for now I'm going to just O(N^{2}) brute force it.
        print "Determining neighbours..."
        for chromophore1ID in self.chromophores:
            for chromophore2ID in self.chromophores:
                if chromophore1ID == chromophore2ID:
                    continue
                if helperFunctions.calculateSeparation(self.chromophores[chromophore1ID]['COMPosn'], self.chromophores[chromophore2ID]['COMPosn']) <= self.maximumHoppingDistance:
                    if self.chromophores[chromophore1ID]['chromoID'] not in self.chromophores[chromophore2ID]['neighbours']:
                        self.chromophores[chromophore2ID]['neighbours'].append(self.chromophores[chromophore1ID]['chromoID'])
                    if self.chromophores[chromophore2ID]['chromoID'] not in self.chromophores[chromophore1ID]['neighbours']:
                        self.chromophores[chromophore1ID]['neighbours'].append(self.chromophores[chromophore2ID]['chromoID'])


    def obtainHoppingPairs(self):
        chromophorePairs = []
        chromophorePairIDs = []
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
        print "There are", len(chromophorePairIDs), "pairs of chromophores"
        for chromophorePair in chromophorePairIDs:
            chromophorePairs.append([self.chromophores[chromophorePair[0]], self.chromophores[chromophorePair[1]]])
        return chromophorePairs







def plotMolecule3D(molecule, moleculeBackbone):
    allAtoms = np.array(molecule['position'])
    COMs = []
    normals = []
    planes = []
    for thioRing in moleculeBackbone:
        normals.append([])
        planes.append([])
        COMs.append(np.array(thioRing['COM']))
        for i in range(10):
            normals[-1].append(np.array(thioRing['COM']+(np.array(thioRing['normal'])*i/5.)))
            planes[-1].append(np.array(thioRing['COM']+(np.array(thioRing['plane'])*i/5.)))
    COMs = np.array(COMs)
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    #ax.scatter(allAtoms[:,0], allAtoms[:,1], allAtoms[:,2], s = 20, c = 'g')
    ax.scatter(COMs[:,0], COMs[:,1], COMs[:,2], s = 50, c = 'r')
    for thioRing in normals:
        for coords in thioRing:
            ax.scatter(coords[0], coords[1], coords[2], s = 20, c = 'b')
    for thioPlane in planes:
        for coords in thioPlane:
            ax.scatter(coords[0], coords[1], coords[2], s = 20, c = 'g')
    #ax.set_xlim((-8,8))
    ax.set_xlim((-40,40))
    ax.set_ylim((-40,40))
    #ax.set_zlim((-15,15))
    ax.set_zlim((-40,40))
    plt.show()






    

def generateORCAInput(self):
    # CreateName (rounded to .1 Ang)
    self.createName()
    # Check that file with the right name doesn't already exist
    #       If it does, pass.
    #       Otherwise make the ORCA inputFile
    exists = self.makeDirs()
    if exists == True:
        print "File", self.fullPath, "already exists, skipping...\n"
        return
    atomsToWrite = self.getAtomsToWrite()
    self.writeInpFile(atomsToWrite)
    # RUN ORCA
    # Analyse ORCA outputs to create a structure (dictionary?) that takes two segment numbers and returns the transferIntegral        



def writeInpFile(self, atomsToWrite):
    templateFile = open('./templates/template.inp', 'r')
    templateLines = templateFile.readlines()
    templateFile.close()
    linesToWrite = []
    for lineNo in range(len(templateLines)):
        if templateLines[lineNo] == '*':
            for atom in atomsToWrite:
                lineToWrite = ' '
                lineToWrite += str(atom[0])+' '
                lineToWrite += str(atom[1])+' '
                lineToWrite += str(atom[2])+' '
                lineToWrite += str(atom[3])+'\n'
                linesToWrite.append(lineToWrite)
        linesToWrite.append(templateLines[lineNo])
    orcaInputFile = open(self.fullPath, 'w+')
    orcaInputFile.writelines(linesToWrite)
    orcaInputFile.close()
    print "Orca Input file written to", self.fullPath, "\n"
