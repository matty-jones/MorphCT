import numpy as np
import helperFunctions
import copy



class obtain:
    def __init__(self, morphologyData, moleculeProperties, boxSize):
        self.morphologyData = morphologyData
        self.moleculeProperties = moleculeProperties
        self.boxSize = boxSize
        self.simDims = [[-boxSize[0]/2.0, boxSize[0]/2.0], [-boxSize[1]/2.0, boxSize[1]/2.0], [-boxSize[2]/2.0, boxSize[2]/2.0]]
        print self.simDims
        self.maximumHoppingDistance = 10.0
        self.getChromoPosns()
        self.updateNeighbourList()
        self.obtainHoppingPairs()


    def getChromoPosns(self):
        self.chromophores = {}
        globalChromoNo = -1 # This number includes any periodic chromophores that are not in the original simulation volume
        self.numberOfPeriodicChromophores = 0
        for molNo, molecule in enumerate(self.moleculeProperties):
            for chromoNo, chromophore in enumerate(molecule['morphologyChromophores']):
                globalChromoNo += 1
                chromoDict = {'molID': molNo, 'chromoID': globalChromoNo, 'periodic': False, 'realChromoID': globalChromoNo, 'position': [], 'image': 0, 'type': [], 'mass': [], 'COMPosn': 0, 'neighbours': []}
                for atomID in chromophore:
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
        for chromoID in self.chromophores:
            for neighbour in self.chromophores[chromoID]['neighbours']:
                if self.chromophores[neighbour]['realChromoID'] == self.chromophores[chromoID]['realChromoID']:
                    print "Chromophore 1 =", self.chromophores[chromoID]
                    print "Chromophore 2 =", self.chromophores[neighbour]
                    print "ChromoIDs =", self.chromophores[chromoID]['chromoID'], self.chromophores[neighbour]['chromoID']
                    print "Neighbours =", self.chromophores[chromoID]['neighbours']
                    raise SystemError("Hopping from one chromo to its own periodic image")
                elif neighbour < chromoID:
                    if [chromoID, neighbour] not in chromophorePairs:
                        chromophorePairs.append([chromoID, neighbour])
                else:
                    if [neighbour, chromoID] not in chromophorePairs:
                        chromophorePairs.append([neighbour, chromoID])
        print "There are", len(chromophorePairs), "pairs of chromophores"

        






def findSegmentPairs(segmentsMaster, atomMaster, maximumHoppingDistance, simVolData):
    # len(segmentsMaster) = number of segments
    print simVolData
    periodicMaster = findPeriodicSegmentsList(segmentsMaster, maximumHoppingDistance, simVolData) # The same structure as segmentsMaster, but only including segments that are not in the initial bounded volume (i.e. are periodically linked to the main system)
    # These periodic segments have the same segment number as in the main, non-image volume, so as long as we use that name for the orca input files, but use the coordinates given in periodicMaster, then we should be ok to calculate transfer integrals.
    # Because of this, treat `periodic pairs' separately to the main segment pairs to ensure we don't get mixed up.
    periodicMasterDict = {} # Useful to sort the list by the segment number for the periodic segments because it no longer lines up with the list indices
    periodicAtomMaster = {}
    for periodicSegment in periodicMaster:
        for atom in periodicSegment.iteritems():
            periodicAtomMaster[atom[0]] = atom[1]

    segmentPairs = []
    periodicSegmentPairs = []
    # Just look at first one first
    for i in range(len(segmentsMaster)):
        segment = segmentsMaster[i]
        nearbySegments = []
        nearbyPeriodicSegments = []
        for currentAtom in segment.iteritems():
            currentAtomCoords = [currentAtom[1][0], currentAtom[1][1], currentAtom[1][2]]
            currentSegment = currentAtom[1][4]
            for compareAtom in atomMaster.iteritems():
                if (compareAtom[1][4] in nearbySegments) or (compareAtom[1][4] == currentSegment) or (currentAtom[0] == compareAtom[0]):
                    continue
                compareAtomCoords = [compareAtom[1][0], compareAtom[1][1], compareAtom[1][2]]
                separation = calcSeparation(currentAtomCoords, compareAtomCoords)
                if separation >= maximumHoppingDistance: # Arbitrarily discard atoms more than maximumHoppingDistance angstroms away
                    continue
                # If it's not continued already:
                # Atom's home segment is not classed as nearby to currentSegment and is within 10 angstroms
                nearbySegments.append(compareAtom[1][4])
            for compareAtom in periodicAtomMaster.iteritems():
                if (compareAtom[1][4] in nearbySegments) or (compareAtom[1][4] == currentSegment) or (currentAtom[0] == compareAtom[0]):
                    continue
                compareAtomCoords = [compareAtom[1][0], compareAtom[1][1], compareAtom[1][2]]
                separation = calcSeparation(currentAtomCoords, compareAtomCoords)
                if separation >= maximumHoppingDistance:
                    continue
                if compareAtom[1][4] not in nearbyPeriodicSegments:
                    nearbyPeriodicSegments.append(compareAtom[1][4])

        if [currentSegment] not in segmentPairs:
            segmentPairs.append([currentSegment])

        # Remove the current segment from the nearby list
        for i in range(len(nearbySegments)):
            if nearbySegments[i] == currentSegment:
                continue
            if currentSegment < nearbySegments[i]:
                segmentPair = [currentSegment, nearbySegments[i]]
            else:
                segmentPair = [nearbySegments[i], currentSegment]
            if segmentPair not in segmentPairs:
                segmentPairs.append(segmentPair)

        for i in range(len(nearbyPeriodicSegments)):
            if nearbyPeriodicSegments[i] == currentSegment:
                continue
            periodicSegmentPair = [currentSegment, nearbyPeriodicSegments[i]] # Always real segment first, periodic second
            if periodicSegmentPair not in periodicSegmentPairs:
                periodicSegmentPairs.append(periodicSegmentPair)

        print "Current Segment =", currentSegment, "Nearby `real' Segments =", nearbySegments, "Nearby `periodic' Segments =", nearbyPeriodicSegments

    # print "SegmentPairs =", segmentPairs
    print "\nLen SegmentPairs =", len(segmentPairs)
    print "Len periodicSegmentPairs =", len(periodicSegmentPairs)
    print "Job done."

    for segment in periodicMaster:
        segmentNumber = segment.items()[0][1][4] # The number for this segment
        periodicMasterDict[segmentNumber] = segment

    return segmentPairs, periodicSegmentPairs, periodicMasterDict



# -------------==================== ======================-----------------




# -------------==================== Creating ORCA Inputs  ======================-----------------

class ORCAInput:
    def __init__(self, segmentsMaster, segmentPair, datName, rotateThio, rotateAlk1, rotateAlk2, justPlotThios, periodicSegmentPair = False, periodicMasterDict = None):
        self.rotateThio = rotateThio
        self.rotateAlk1 = rotateAlk1
        self.rotateAlk2 = rotateAlk2
        self.justPlotThios = justPlotThios
        self.segmentPair = segmentPair
        self.datName = datName[:-4]
        self.segmentAtoms = []
        if periodicSegmentPair == False:
            self.segmentAtoms.append(segmentsMaster[segmentPair[0]-1])
            COMSeg1 = self.determineCOMFromCG(segmentsMaster[segmentPair[0]-1])
            COMSeg2 = None
            if len(segmentPair) == 2:
                COMSeg2 = self.determineCOMFromCG(segmentsMaster[segmentPair[1]-1])
                self.segmentAtoms.append(segmentsMaster[segmentPair[1]-1])
            self.writeCOMCSV(COMSeg1, COMSeg2)
        else:
            self.segmentAtoms.append(segmentsMaster[segmentPair[0]-1])
            self.segmentAtoms.append(periodicMasterDict[segmentPair[1]])
            COMSeg1 = self.determineCOMFromCG(segmentsMaster[segmentPair[0]-1])
            COMSeg2 = self.determineCOMFromCG(periodicMasterDict[segmentPair[1]])
            self.writeCOMCSV(COMSeg1, COMSeg2)
        self.monomerAtoms = self.read3HTTemplate() # These are the coordinates of the atoms representing a CG site at [0, 0, 0]

    def writeCOMCSV(self, COMData1, COMData2):
        # Read the current data into a dictionary to make sure that we're not re-writing the same data
        currentSegData = {}
        if COMData2 == None:
            COMCSVName = './orcaOutputs/'+self.datName+'/SegmentSingleCOM.csv'
            csvReadOnly = open(COMCSVName, 'a+')
            csvReadOnly.seek(0)
            csvData = csv.reader(csvReadOnly, delimiter = ',')
            for row in csvData:
                currentSegData[row[0]] = 1
            csvReadOnly.close()
            currentSegment = COMData1[0]
            if (currentSegment not in currentSegData):
                csvFile = open(COMCSVName, 'a')
                csvWriter = csv.writer(csvFile, delimiter = ',')
                csvWriter.writerow([COMData1[0], COMData1[1][0], COMData1[1][1], COMData1[1][2]])
                csvFile.close()
        else:
            COMCSVName = './orcaOutputs/'+self.datName+'/SegmentPairCOM.csv'
            csvReadOnly = open(COMCSVName, 'a+')
            csvReadOnly.seek(0)
            csvData = csv.reader(csvReadOnly, delimiter = ',')
            for row in csvData:
                currentSegData[(row[0], row[1])] = 1
            csvReadOnly.close()
            currentSegmentPair = (COMData1[0], COMData2[0])
            if (currentSegmentPair not in currentSegData):
                separation = calcSeparation([COMData1[1][0], COMData1[1][1], COMData1[1][2]], [COMData2[1][0], COMData2[1][1], COMData2[1][2]])
                csvFile = open(COMCSVName, 'a')
                csvWriter = csv.writer(csvFile, delimiter = ',')
                # Write the row as [seg1ID, seg2ID, seg1X, seg1Y, seg1Z, seg2X, seg2Y, seg2Z, separation]
                csvWriter.writerow([COMData1[0], COMData2[0], COMData1[1][0], COMData1[1][1], COMData1[1][2], COMData2[1][0], COMData2[1][1], COMData2[1][2], separation])
                csvFile.close()

    def determineCOMFromCG(self, inputAtoms):
        # SegmentsMaster-like object comes in
        massOfThio = 81.11657
        massOfAlk1 = 42.08127
        massOfAlk2 = 43.08924
        massWeightedX = 0.
        massWeightedY = 0.
        massWeightedZ = 0.
        totalMass = 0.
        for atomData in inputAtoms.values():
            currentSegment = atomData[4]

            massWeightedX += atomData[0]*massOfThio
            massWeightedY += atomData[1]*massOfThio
            massWeightedZ += atomData[2]*massOfThio
            totalMass += massOfThio

            massWeightedX += atomData[8][0]*massOfAlk1
            massWeightedY += atomData[8][1]*massOfAlk1
            massWeightedZ += atomData[8][2]*massOfAlk1
            totalMass += massOfAlk1

            massWeightedX += atomData[9][0]*massOfAlk2
            massWeightedY += atomData[9][1]*massOfAlk2
            massWeightedZ += atomData[9][2]*massOfAlk2
            totalMass += massOfAlk2
        return [currentSegment, np.array([massWeightedX/float(totalMass), massWeightedY/float(totalMass), massWeightedZ/float(totalMass)])]


    def determineCOM(self, inputAtoms, thioOnly = False, alk1Only = False, alk2Only = False):
        massWeightedX = 0.
        massWeightedY = 0.
        massWeightedZ = 0.
        totalMass = 0.
        atoms = []
        if thioOnly == True: # These are the lines that contain thiophene atoms in the .xyz
            atoms.append(inputAtoms[0]) # S1
            atoms.append(inputAtoms[1]) # C2
            atoms.append(inputAtoms[3]) # C4
            atoms.append(inputAtoms[4]) # C5
            atoms.append(inputAtoms[5]) # C6
            atoms.append(inputAtoms[7]) # H8
        elif alk1Only == True: # These are the lines that contain alk1 atoms in the .xyz
            atoms.append(inputAtoms[8]) # C9
            atoms.append(inputAtoms[10]) # C11
            atoms.append(inputAtoms[13]) # C14
            atoms.append(inputAtoms[12]) # H13
            atoms.append(inputAtoms[2]) # H3
            atoms.append(inputAtoms[7]) # H6
            atoms.append(inputAtoms[15]) # H16
            atoms.append(inputAtoms[9]) # H10
            atoms.append(inputAtoms[17]) # H18
        elif alk2Only == True: # These are the lines that contain alk2 atoms in the .xyz
            atoms.append(inputAtoms[16])
            atoms.append(inputAtoms[11])
            atoms.append(inputAtoms[20])
            atoms.append(inputAtoms[19])
            atoms.append(inputAtoms[14])
            atoms.append(inputAtoms[22])
            atoms.append(inputAtoms[21])
            atoms.append(inputAtoms[18])
            atoms.append(inputAtoms[23])
            atoms.append(inputAtoms[24])
        else:
            atoms = inputAtoms
        for atom in atoms:
            if atom[0] == 0: # Determining chain COM, therefore mass factor is irrelevent
                mass = 1.0                
            elif atom[0] == 'S':
                mass = 32.065
            elif atom[0] == 'C':
                mass = 12.0107
            elif atom[0] == 'H':
                mass = 1.00794

            massWeightedX += atom[1][0]*mass
            massWeightedY += atom[1][1]*mass
            massWeightedZ += atom[1][2]*mass
            totalMass += mass
        return np.array([massWeightedX/float(totalMass), massWeightedY/float(totalMass), massWeightedZ/float(totalMass)])

    def read3HTTemplate(self):
        monomerFile = open('./templates/3HT_topdown.xyz', 'r')
        monomer = monomerFile.readlines()
        monomerFile.close()
        # Check centrepoint
        atoms = []
        for line in monomer[2:-1]:
            tempLine = []
            coordinates = []
            for character in line.split(' '):
                if len(character) != 0:
                    if len(character) == 1: # This is the element of the atom
                        tempLine.append(str(character))
                    else: # This is a coordinate
                        try:
                            coordinates.append(float(character))
                        except: # "\n" is breaking it
                            coordinates.append(float(character[:-2]))
            tempLine.append(coordinates)
            atoms.append(tempLine)

        COM = self.determineCOM(atoms, thioOnly = True)
        for atom in atoms:
            coordArray = np.array(atom[1])
            coordArray -= COM
            atom[1] = list(coordArray)

        return atoms


    def createName(self):
        # Name was too complicated, made it easier
        if (len(segmentPair) == 1):
            segNumber = str(segmentPair[0])
            while len(segNumber) < 4:
                segNumber = '0'+segNumber
            self.name = 'seg'+segNumber+'.inp'
        elif (len(segmentPair) == 2):
            segNumber1 = str(segmentPair[0])
            while len(segNumber1) < 4:
                segNumber1 = '0'+segNumber1
            segNumber2 = str(segmentPair[1])
            while len(segNumber2) < 4:
                segNumber2 = '0'+segNumber2
            self.name = 'seg'+segNumber1+'_seg'+segNumber2+'.inp'
        else:
            print "Segment Pair =", segmentPair
            raise SystemError("Segment Pair length is either 0 or > 2")


        # # Want a name kinda like C1_LengthOfChain1_C2_LengthOfChain2_X_XSeparationFromCentrePoint_Y_YSep_Z_ZSep.inp
        # C1Length = len(self.segmentAtoms[0])
        # C2Length = len(self.segmentAtoms[1])
        # # Determine COMs of each chain
        # chain1 = []
        # chain2 = []
        # for atom in self.segmentAtoms[0].iteritems():
        #     chain1.append([0, [atom[1][0], atom[1][1], atom[1][2]]])
        # for atom in self.segmentAtoms[1].iteritems():
        #     chain2.append([0, [atom[1][0], atom[1][1], atom[1][2]]])
        # COM1 = self.determineCOM(chain1)
        # COM2 = self.determineCOM(chain2)
        # separation = COM2-COM1
        # self.name = 'C1_'+str(len(chain1))+'_C2_'+str(len(chain2))+'_X_%.1f_Y_%.1f_Z_%.1f' % (float(separation[0]), float(separation[1]), float(separation[2]))


    def makeDirs(self):
        # Make the correct directory
        dirs = os.listdir('./orcaInputs')
        if self.datName not in dirs:
            os.makedirs('./orcaInputs/'+self.datName)
        filesList = os.listdir('./orcaInputs/'+self.datName)
        self.fullPath = './orcaInputs/'+self.datName+'/'+self.name
        if self.name in filesList:
            return True
        else:
            return False

    def determineThioToAlk1Vector(self, atoms):
        #thioPosn = [0, 0, 0] because we moved the monomer to the middle when we read it
        thioPosn = self.determineCOM(atoms, thioOnly = True)
        alk1Posn = self.determineCOM(atoms, alk1Only = True)
        # print "ThioPosn =", thioPosn
        # print "alk1Posn =", alk1Posn
        alk1Vector = normaliseVec(alk1Posn-thioPosn)
        return alk1Vector



    def determineNormalToThiopheneRing(self, monomerAtoms):
        # Plane is determined by two vectors, one from the Sulphur (S1 - atom 1 in the .xyz) to a neighbouring Carbon (C4 - atom 4 in the .xyz),
        # and the other between the sulphur (S1 - atom 1) to a carbon on the opposite side of the ring (C5 - atom 5 in the .xyz)
#        sulphurAtom = monomerAtoms[0]

        firstCarbon = monomerAtoms[1]
        neighbourCarbon = monomerAtoms[4]
        oppositeCarbon = monomerAtoms[3]
        
        vec1 = findAxis(firstCarbon[1], neighbourCarbon[1])
        vec2 = findAxis(firstCarbon[1], oppositeCarbon[1])


        normalToRing = np.cross(vec2, vec1)
        # normalToRingMagnitude = np.sqrt(((normalToRing[0])**2) + ((normalToRing[1])**2) + ((normalToRing[2])**2))
        # normalToRing /= normalToRingMagnitude
        normalToRing = normaliseVec(normalToRing)

        return normalToRing


    # def determineNormalToThiopheneRing(self, monomerAtoms, thioAlk1Vector):
    #     # Plane is determined by two vectors, one from the thiophene-alk1 vector (given)
    #     # And the other between the atoms C2 and C4 (which are the carbons that bond to the adjacent monomers)
    #     print monomerAtoms[1], monomerAtoms[3]
    #     leftCarbon = monomerAtoms[1][1] # C2
    #     rightCarbon = monomerAtoms[3][1] # C4

    #     crossRingVector = findAxis(leftCarbon, rightCarbon)

    #     print crossRingVector, "\n"

    #     normalToRing = np.cross(crossRingVector, thioAlk1Vector)
    #     normalToRing = normaliseVec(normalToRing)

    #     return normalToRing

    def rotationMatrixAroundAxis(self, theta, axis):
        # Rotation matrix obtained from http://goo.gl/RkW80
        R = np.matrix([[ np.cos(theta)+(axis[0]**2)*(1-np.cos(theta)), axis[0]*axis[1]*(1-np.cos(theta))-axis[2]*np.sin(theta), axis[0]*axis[2]*(1-np.cos(theta))+axis[1]*np.sin(theta) ],
                       [ axis[1]*axis[0]*(1-np.cos(theta))+axis[2]*np.sin(theta), np.cos(theta)+(axis[1]**2)*(1-np.cos(theta)), axis[1]*axis[2]*(1-np.cos(theta))-axis[0]*np.sin(theta) ],
                       [ axis[2]*axis[0]*(1-np.cos(theta))-axis[1]*np.sin(theta), axis[2]*axis[1]*(1-np.cos(theta))+axis[0]*np.sin(theta), np.cos(theta)+(axis[2]**2)*(1-np.cos(theta)) ]])
        return R


    def executeRotation(self, testRotationAngle, rotationAxis, leftCarbon, rightCarbon):
        leftCarbon = copy.deepcopy(leftCarbon)
        rightCarbon = copy.deepcopy(rightCarbon)
        rotationMatrix = self.rotationMatrixAroundAxis(testRotationAngle, rotationAxis)
        newLeftCarbonCoords = np.array(np.transpose(rotationMatrix*np.transpose(np.matrix(leftCarbon))))[0].tolist()
        newRightCarbonCoords = np.array(np.transpose(rotationMatrix*np.transpose(np.matrix(rightCarbon))))[0].tolist()
        return newLeftCarbonCoords, newRightCarbonCoords


    def calcSeparationAngle(self, vec1, vec2):
        vec1 = normaliseVec(vec1)
        vec2 = normaliseVec(vec2)
        dotProduct = np.dot(vec1, vec2)
        theta = np.arccos(dotProduct)
        theta = abs(theta)
        if theta > np.pi:
            theta = theta-np.pi
        return theta



    def thioInPlaneRotation(self, atoms, adjThioCoords, rotationAxis):
        # print "\n"
        totalAngleTurnedBy = 0
        testRotationAngle = np.pi/180.
        # Target vector is the vector between the two neighbour thiophenes, or just the vector to the adjacent thio if only one bonded neighbour
        if len(adjThioCoords) == 2:
            targetVector = findAxis(adjThioCoords[0], adjThioCoords[1])
        else:
            # Ends of chains don't really work because they only have one neighbour. As such, let's ignore them entirely and do not rotate
            return np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        # Current vector is the vector between the left and right carbons in the chain
        leftCarbon = copy.deepcopy(atoms[3][1])
        rightCarbon = copy.deepcopy(atoms[1][1])
        # print "Left Carbon", leftCarbon, "Right Carbon", rightCarbon
        currentVector = findAxis(leftCarbon, rightCarbon)
        currentTheta = self.calcSeparationAngle(currentVector, targetVector)
        # To minimise, perform rotation. If arccos of dot product > pi, subtract pi. (just means going opposite way so that's fine). Minimise theta.
        # First, perform a test rotation to see if we're rotating the right way
        # print "Target vector =", targetVector
        # print "Current vector =", currentVector
        # print "Beginning test rotation..."
        newLeftCarbonCoords, newRightCarbonCoords = self.executeRotation(testRotationAngle, rotationAxis, leftCarbon, rightCarbon)
        newVector = findAxis(newLeftCarbonCoords, newRightCarbonCoords)
        newTheta = self.calcSeparationAngle(newVector, targetVector)
        # print "Current Theta =", currentTheta, "new Theta =", newTheta
        if newTheta > currentTheta:
            # print "Rotating the wrong way! Reversing test rotation..."
            # Rotating the wrong way!
            testRotationAngle = -testRotationAngle
            newLeftCarbonCoords, newRightCarbonCoords = self.executeRotation(testRotationAngle, rotationAxis, leftCarbon, rightCarbon)
            newVector = findAxis(newLeftCarbonCoords, newRightCarbonCoords)
            newTheta = self.calcSeparationAngle(newVector, targetVector)
            if newTheta > currentTheta:
                print "Got an issue, rotated the opposite way but it didn't help."
                print "Likely means that we are already at the optimal position"
                return np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        # newTheta is now < currentTheta
        # print "---=== NEW VECTOR =", newVector, "===---"
        while newTheta <= currentTheta:
            # print "Beginning of loop: New Theta =", newTheta, "Current Theta =", currentTheta
            currentTheta = copy.deepcopy(newTheta)
            currentVector = copy.deepcopy(newVector)
            totalAngleTurnedBy += testRotationAngle
            # print "Executing another rotation (currentTheta =", str(currentTheta)+", totalAngleTurnedBy =", str(totalAngleTurnedBy)+")"
            newLeftCarbonCoords, newRightCarbonCoords = self.executeRotation(testRotationAngle, rotationAxis, leftCarbon, rightCarbon)
            # print "Left Carbon =", leftCarbon, "newLeftCarbon =", newLeftCarbonCoords
            # print "Right Carbon =", rightCarbon, "newRightCarbon =", newRightCarbonCoords
            newVector = findAxis(newLeftCarbonCoords, newRightCarbonCoords)
            newTheta = self.calcSeparationAngle(newVector, targetVector)
            # print "NEW VECTOR =", newVector
            # print "End of loop: Current theta is still", currentTheta, " but now newTheta =", newTheta, "(which has to be < currentTheta)"
            leftCarbon = newLeftCarbonCoords
            rightCarbon = newRightCarbonCoords

        rotMatrix = self.rotationMatrixAroundAxis(totalAngleTurnedBy, rotationAxis)


        # print "NewLeftCarbon", newLeftCarbonCoords, "NewRightCarbon", newRightCarbonCoords
        # print "NewVector =", newVector

        return rotMatrix
                


    def rotationMatrix(self, vector1, vector2):
        # A function to return the rotation matrix around the origin that maps vector1 to vector 2
        crossProduct = np.cross(vector1, vector2)
        sinAngle = np.sqrt(((crossProduct[0]**2) + ((crossProduct[1])**2) + ((crossProduct[2])**2)))
        cosAngle = np.dot(vector1, vector2)

        skewMatrix = np.matrix([[0, -crossProduct[2], crossProduct[1]], [crossProduct[2], 0, -crossProduct[0]], [-crossProduct[1], crossProduct[0], 0]])
        skewMatrixSquared = skewMatrix * skewMatrix
        
        rotMatrix = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) + skewMatrix + skewMatrixSquared*((1 - cosAngle)/(sinAngle**2))

        return rotMatrix



    def manipulateMonomer(self, COMThio, CGThioPlaneNormal, CGAdjThioCoords, CGThioAlk1, COMAlk1, COMAlk2):
        alk1ManipulatedMonomer = []
        alk1AndThioManipulatedMonomer = []
        alk1AndThioZManipulatedMonomer = []
        thioAlk1CompletedMonomer = []
        thioAlk1Alk2CompletedMonomer = []
        monomerFinal = []

        relativeAlk1Posn = [(COMAlk1[0] - COMThio[0]), (COMAlk1[1] - COMThio[1]), (COMAlk1[2] - COMThio[2])] # Alk1 posn with thio as origin
        relativeAlk2Posn = [(COMAlk2[0] - COMThio[0]), (COMAlk2[1] - COMThio[1]), (COMAlk2[2] - COMThio[2])] # Alk1 posn with thio as origin

        templateAlk1Vector = self.determineThioToAlk1Vector(self.monomerAtoms)

        R_alk1 = self.rotationMatrix(templateAlk1Vector, CGThioAlk1)

        for atom in self.monomerAtoms:
            transformedAtom = copy.deepcopy(atom)
            transformedAtomCoords = np.transpose(R_alk1*np.transpose(np.matrix(transformedAtom[1])))
            transformedAtom[1] = [transformedAtomCoords[0,0], transformedAtomCoords[0,1], transformedAtomCoords[0,2]]
            alk1ManipulatedMonomer.append(transformedAtom)


        thioNormalVector = self.determineNormalToThiopheneRing(alk1ManipulatedMonomer)
#        thioNormalVector = self.determineNormalToThiopheneRing(alk1ManipulatedMonomer, CGThioAlk1) # We just mapped the templateAlk1Vector to CGThioAlk1, so its Alk1 vector should already be CGThioAlk1

        R_thio = self.rotationMatrix(thioNormalVector, CGThioPlaneNormal)

        # Try only rotating the thiophene ring and leave the alkyl sidechain in its place
        thioRingElements = [0, 1, 3, 4, 5, 7] # These are the line numbers of the atoms that belong to the thiophene ring
        alk1Elements = [8, 2, 12, 10, 6, 15, 13, 9, 17] # These line numbers are the alk1 atoms
        alk2Elements = [16, 11, 20, 19, 14, 22, 21, 18, 23, 24] # These line numbers are the alk2 atoms

        for atomNo in range(len(alk1ManipulatedMonomer)):
            transformedAtom = copy.deepcopy(alk1ManipulatedMonomer[atomNo])
            transformedAtomCoords = np.transpose(R_thio*np.transpose(np.matrix(transformedAtom[1])))
            transformedAtom[1] = [transformedAtomCoords[0,0], transformedAtomCoords[0,1], transformedAtomCoords[0,2]]
            alk1AndThioManipulatedMonomer.append(transformedAtom)

        # # ORCA inputs still don't work, so now we need to rotate the ring around the normal to this ring

        newThioNormalVector = self.determineNormalToThiopheneRing(alk1AndThioManipulatedMonomer)
        
        # for atom in alk1AndThioManipulatedMonomer:
        #     alk1AndThioFullyManipulatedMonomer.append(atom)


        if self.rotateThio == True:

            R_thioZ = self.thioInPlaneRotation(alk1AndThioManipulatedMonomer, CGAdjThioCoords, newThioNormalVector)

            for atomNo in range(len(alk1AndThioManipulatedMonomer)):
                transformedAtom = copy.deepcopy(alk1AndThioManipulatedMonomer[atomNo])
                if atomNo in thioRingElements: #  Just Rotate the Thiophene Ring
                    transformedAtomCoords = np.transpose(R_thioZ*np.transpose(np.matrix(transformedAtom[1])))
                    transformedAtom[1] = [transformedAtomCoords[0,0], transformedAtomCoords[0,1], transformedAtomCoords[0,2]]
                alk1AndThioZManipulatedMonomer.append(transformedAtom)
        else:
            for atom in alk1AndThioManipulatedMonomer:
                alk1AndThioZManipulatedMonomer.append(atom)



        if self.rotateAlk1 == True:

            # ---=== How to fix the rest of the geometry: ===---
            # 0) Calculate the current COM position of the Alk1 chain (C9, H3, H13, C11, H7, H16, C14, H10, H18)
            currentAlk1COM = self.determineCOM(alk1AndThioZManipulatedMonomer, alk1Only = True) # As if the origin is the thio
            # 0.5) Move all of the Alk1 atoms so that their new COM position matches that of COMAlk1
            alk1Movement = [(relativeAlk1Posn[0] - currentAlk1COM[0]), (relativeAlk1Posn[1] - currentAlk1COM[1]), (relativeAlk1Posn[2] - currentAlk1COM[2])]
            # 1) Make a coordinate switch so that COMAlk1 is now the origin (Move all atoms st. CurrentCoordinates = CurrentCoordinate - relativeAlk1Posn)
            for atomNo in range(len(alk1AndThioZManipulatedMonomer)):
                if atomNo in alk1Elements:
                    alk1AndThioZManipulatedMonomer[atomNo][1] = list(np.array(alk1AndThioZManipulatedMonomer[atomNo][1]) + np.array(alk1Movement))
                alk1AndThioZManipulatedMonomer[atomNo][1] = list(np.array(alk1AndThioZManipulatedMonomer[atomNo][1]) - np.array(relativeAlk1Posn))
            # 2) "CurrentVector" = findAxis( C14, C9 )
            currentAlk1Vector = findAxis(alk1AndThioZManipulatedMonomer[13][1], alk1AndThioZManipulatedMonomer[8][1])
            # 3) "TargetVector" = findAxis( COMAlk1, C6 ) # NB COMAlk1 == 0
            # Might not be C6, sometimes ThioZ flips the pentagon but as it is regular, we can just find the closest carbon out of C5, C6
            C5Sep = calcSeparation(relativeAlk1Posn, alk1AndThioZManipulatedMonomer[4][1])
            C6Sep = calcSeparation(relativeAlk1Posn, alk1AndThioZManipulatedMonomer[5][1])
            if (C5Sep < C6Sep):
                # Connect to C5
                targetVector = findAxis(np.array([0,0,0]), alk1AndThioZManipulatedMonomer[4][1])
                # targetVector = findAxis(np.array(relativeAlk1Posn), alk1AndThioZManipulatedMonomer[4][1])
            else:
                # Connect to C6
                targetVector = findAxis(np.array([0,0,0]), alk1AndThioZManipulatedMonomer[5][1])
                # targetVector = findAxis(np.array(relativeAlk1Posn), alk1AndThioZManipulatedMonomer[5][1])
            # 4) Map currentVector onto targetVector as before
            R_alk1_rot = self.rotationMatrix(currentAlk1Vector, targetVector)

            for atomNo in range(len(alk1AndThioZManipulatedMonomer)):
                transformedAtom = copy.deepcopy(alk1AndThioZManipulatedMonomer[atomNo])
                if atomNo in alk1Elements:
                    transformedAtomCoords = np.transpose(R_alk1_rot*np.transpose(np.matrix(transformedAtom[1])))
                    transformedAtom[1] = [transformedAtomCoords[0,0], transformedAtomCoords[0,1], transformedAtomCoords[0,2]]
                # 5) Move all atoms by relativeAlk1Posn (current += rel) to return the ThioCOM to the origin again
                transformedAtom[1] = list(np.array(transformedAtom[1]) + np.array(relativeAlk1Posn))
                thioAlk1CompletedMonomer.append(transformedAtom)

        else:
            for atom in alk1AndThioZManipulatedMonomer:
                thioAlk1CompletedMonomer.append(atom)



        if self.rotateAlk2 == True:


            # A) REPEAT ALL STEPS AGAIN FOR THE ALK2 CHAIN
            # 0) Calculate the current COM position of the Alk2 chain = ( C17, H12, H21, C20, H15, H23, C22, H19, H24, H25 )
            currentAlk2COM = self.determineCOM(thioAlk1CompletedMonomer, alk2Only = True) # As if the origin is the thio
            # 0.5) Move all of the Alk2 atoms so that their new COM position matches that of COMAlk1
            alk2Movement = [(relativeAlk2Posn[0] - currentAlk2COM[0]), (relativeAlk2Posn[1] - currentAlk2COM[1]), (relativeAlk2Posn[2] - currentAlk2COM[2])]
            # 1) Make a coordinate switch so that COMAlk2 is now the origin (Move all atoms st. CurrentCoordinates = CurrentCoordinate - relativeAlk2Posn)
            for atomNo in range(len(thioAlk1CompletedMonomer)):
                if atomNo in alk2Elements:
                    thioAlk1CompletedMonomer[atomNo][1] = list(np.array(thioAlk1CompletedMonomer[atomNo][1]) + np.array(alk2Movement))
                thioAlk1CompletedMonomer[atomNo][1] = list(np.array(thioAlk1CompletedMonomer[atomNo][1]) - np.array(relativeAlk2Posn))
            # 2) "CurrentVector" = findAxis( C22, C17 )
            currentAlk2Vector = findAxis(thioAlk1CompletedMonomer[21][1], thioAlk1CompletedMonomer[16][1])
            # 3) "TargetVector" = findAxis( COMAlk2, C14 ) # NB COMAlk2 == 0
            targetVector = findAxis(np.array([0, 0, 0]), thioAlk1CompletedMonomer[13][1])

            # targetVector = findAxis(np.array(relativeAlk2Posn), thioAlk1CompletedMonomer[13][1])
            # 4) Map currentVector onto targetVector as before
            R_alk2_rot = self.rotationMatrix(currentAlk2Vector, targetVector)

            for atomNo in range(len(thioAlk1CompletedMonomer)):
                transformedAtom = copy.deepcopy(thioAlk1CompletedMonomer[atomNo])
                if atomNo in alk2Elements:
                    transformedAtomCoords = np.transpose(R_alk2_rot*np.transpose(np.matrix(transformedAtom[1])))
                    transformedAtom[1] = [transformedAtomCoords[0,0], transformedAtomCoords[0,1], transformedAtomCoords[0,2]]
                # 5) Move all atoms by relativeAlk2Posn (current += rel) to return the ThioCOM to the origin again
                transformedAtom[1] = list(np.array(transformedAtom[1]) + np.array(relativeAlk2Posn))
                thioAlk1Alk2CompletedMonomer.append(transformedAtom)

        else:
            for atom in thioAlk1CompletedMonomer:
                thioAlk1Alk2CompletedMonomer.append(atom)


        # As a final check, somehow after the thioZ rotation, the hydrogen atom ends up being bonded to the wrong carbon and is too close to the alk1 chain.
        # This probably means that the thioZ is rotating the wrong way, but as we're dealing with a regular pentagon, we could quite easily just flip the hydrogen round (as it also makes basically no contribution to the ZINDO measurements, but could possibly prevent the SCF converging if it's located within another atom)
        # To this end, rotate the hydrogen atom (H8) 180 degrees around the axis S1 -> <C6, C5> if the separation between H8 and C9 is < 1A
        separationCheck = calcSeparation(thioAlk1Alk2CompletedMonomer[7][1], thioAlk1Alk2CompletedMonomer[8][1])
        if separationCheck < 1:
            sulphurAtom = thioAlk1Alk2CompletedMonomer[0][1]
            # Opposite side is the average position of C5 and C6
            oppositeSide = [(thioAlk1Alk2CompletedMonomer[4][1][0] + thioAlk1Alk2CompletedMonomer[5][1][0])/2., (thioAlk1Alk2CompletedMonomer[4][1][1] + thioAlk1Alk2CompletedMonomer[5][1][1])/2., (thioAlk1Alk2CompletedMonomer[4][1][2] + thioAlk1Alk2CompletedMonomer[5][1][2])/2.]
            thioSplitAxis = findAxis(sulphurAtom, oppositeSide)
            R_H8_rot = self.rotationMatrixAroundAxis(np.pi, thioSplitAxis)
            transformedAtom = copy.deepcopy(thioAlk1Alk2CompletedMonomer[7])
            transformedAtomCoords = np.transpose(R_H8_rot*np.transpose(np.matrix(transformedAtom[1])))
            transformedAtom[1] = [transformedAtomCoords[0,0], transformedAtomCoords[0,1], transformedAtomCoords[0,2]]
            thioAlk1Alk2CompletedMonomer[7] = transformedAtom


        # Add on an additional Hydrogen Atom
        
        if (self.justPlotThios == True):

            tempHydrogen = copy.deepcopy(thioAlk1Alk2CompletedMonomer[7]) #  No longer atom[7] - it's still H8, but we cut out atoms so it's now just the last one
            sulphurAtom = thioAlk1Alk2CompletedMonomer[0][1]
            oppositeSide = [(thioAlk1Alk2CompletedMonomer[4][1][0] + thioAlk1Alk2CompletedMonomer[5][1][0])/2., (thioAlk1Alk2CompletedMonomer[4][1][1] + thioAlk1Alk2CompletedMonomer[5][1][1])/2., (thioAlk1Alk2CompletedMonomer[4][1][2] + thioAlk1Alk2CompletedMonomer[5][1][2])/2.]
            thioSplitAxis = findAxis(sulphurAtom, oppositeSide)
            R_H8_rot = self.rotationMatrixAroundAxis(np.pi, thioSplitAxis)
            tempHydrogenCoords = np.transpose(R_H8_rot*np.transpose(np.matrix(tempHydrogen[1])))
            tempHydrogen[1] = [tempHydrogenCoords[0,0], tempHydrogenCoords[0,1], tempHydrogenCoords[0,2]]
            thioAlk1Alk2CompletedMonomer.append(tempHydrogen)




        if (self.justPlotThios == True):
            for atomNo in range(len(thioAlk1Alk2CompletedMonomer)):
                if atomNo in thioRingElements: # Just the thiophene ring
                    monomerFinal.append(thioAlk1Alk2CompletedMonomer[atomNo])
            monomerFinal.append(thioAlk1Alk2CompletedMonomer[-1]) # The extra Hydrogen
        else:
            for atom in thioAlk1Alk2CompletedMonomer:
                monomerFinal.append(atom)
        # for atomNo in range(len(thioAlk1Alk2CompletedMonomer)):
        #     # if (atomNo == 7):
        #     #     continue
        #     # if (atomNo in thioRingElements):
        #     testMonomer.append(thioAlk1Alk2CompletedMonomer[atomNo])


        # return testMonomer




        # Finally add on the ThioCOM Coords

        for atom in monomerFinal:
            atom[1] = [atom[1][0]+COMThio[0], atom[1][1]+COMThio[1], atom[1][2]+COMThio[2]]





        return monomerFinal



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
        


    def getAtomsToWrite(self):
        atomsToWrite = []
        segmentNo = 1
        atomNo = 1
        for segment in self.segmentAtoms:
            for CGAtom in segment.iteritems():
                # print "Treating segment", str(segmentNo)+", monomer", str(atomNo)+"."
                CGCoord = np.array([CGAtom[1][0], CGAtom[1][1], CGAtom[1][2]])
                CGThioPlaneNormal = CGAtom[1][5]
                CGThioAlk1Axis = CGAtom[1][6]
                CGAdjThioCoords = CGAtom[1][7]
                CGAlk1 = CGAtom[1][8]
                CGAlk2 = CGAtom[1][9]
                fineGrainedAtoms = self.manipulateMonomer(CGCoord, CGThioPlaneNormal, CGAdjThioCoords, CGThioAlk1Axis, CGAlk1, CGAlk2)
                for atom in fineGrainedAtoms:
                    atomsToWrite.append([atom[0], str(atom[1][0]), str(atom[1][1]), str(atom[1][2])])
                atomNo += 1
            segmentNo += 1
        # Write ORCA File
        return atomsToWrite

