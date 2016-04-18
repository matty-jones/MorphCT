import numpy as np
import copy
import os
import pickle
import multiprocessing as mp
#from cme_utils.manip import pbc

def findMagnitude(vector):
    '''This function simply returns the magnitude of a given vector'''
    return np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)


def updateXMLBoxLength(adjustedInputFileName, boxSize):
    '''This function opens a hoomd xml and updates it to have a given simulation volume (boxSize)'''
    with open(adjustedInputFileName, 'r') as xmlFile:
        xmlData = xmlFile.readlines()
    for lineNo in range(len(xmlData)):
        if 'box' in xmlData[lineNo]:
            newBoxLine = ''
            quoteMarksLoc = findIndex(xmlData[lineNo], '"')
            # The quote marks 0 and 1 are around the number for lx, 2 and 3 are ly,
            # 4 and 5 are lz. Others are for skew (xy, xz, yz)
            listOfLine = list(xmlData[lineNo])
            listOfLine[quoteMarksLoc[4]+1:quoteMarksLoc[5]] = str(boxSize[2])
            listOfLine[quoteMarksLoc[2]+1:quoteMarksLoc[3]] = str(boxSize[1])
            listOfLine[quoteMarksLoc[0]+1:quoteMarksLoc[1]] = str(boxSize[0])
            for character in listOfLine:
                newBoxLine += character
            newBoxLine += '\n'
            xmlData[lineNo] = newBoxLine
            break
    with open(adjustedInputFileName, 'w+') as xmlFile:
        xmlFile.writelines(xmlData)

        
def findIndex(string, character):
    '''This function returns the locations of an inputted character in an inputted string'''
    index = 0
    locations = []
    while index < len(string):
        if string[index] == character:
            locations.append(index)
        index += 1
    if len(locations) == 0:
        return None
    return locations


def calculateSeparation(atom1, atom2):
    '''This function calculates the distance between two input points (either as lists or np.arrays)'''
    xdif = atom1[0] - atom2[0]
    ydif = atom1[1] - atom2[1]
    zdif = atom1[2] - atom2[2]
    return np.sqrt(xdif**2 + ydif**2 + zdif**2)


def linearInterpDescendingY(targetValue, xArray, yArray):
    '''This function takes in two numpy arrays, and then linearly interpolates to find the value of X when Y is equal to targetValue. yArray must be a descending array (doesn't have to be monotonic, but the function will report the first point at which the curve is below targetValue so be careful of noise!). The function returns a value of None if the yArray never drops below the targetValue'''
    xVal = None
    for index, value in enumerate(yArray):
        if value > targetValue:
            continue
        xLo = xArray[index-1]
        xHi = xArray[index]
        yHi = yArray[index-1]
        yLo = yArray[index]
        yDiff = yHi-yLo
        xDiff = xHi-xLo
        yDeltaFrac = (yHi-targetValue)/yDiff
        xVal = xLo + yDeltaFrac*xDiff
        break
    return xVal


def calcCOM(listOfPositions, listOfMasses):
    '''This function calculates the centre of mass of a collection of sites/atoms (listOfPositions) with corresponding mass (listOfMasses)'''
    massWeightedX = 0.
    massWeightedY = 0.
    massWeightedZ = 0.
    totalMass = np.sum(listOfMasses)
    for atomID in range(len(listOfPositions)):
        massWeightedX += listOfPositions[atomID][0]*listOfMasses[atomID]
        massWeightedY += listOfPositions[atomID][1]*listOfMasses[atomID]
        massWeightedZ += listOfPositions[atomID][2]*listOfMasses[atomID]
    return np.array([massWeightedX/float(totalMass), massWeightedY/float(totalMass), massWeightedZ/float(totalMass)])

        
def findAxis(atom1, atom2, normalise=True):
    '''This function determines the normalised vector from the location of atom1 to atom2. The positions can enter as lists or arrays, but are output as arrays'''
    xSep = atom2[0] - atom1[0]
    ySep = atom2[1] - atom1[1]
    zSep = atom2[2] - atom1[2]
    if normalise == True:
        axisVector = normaliseVec(np.array([xSep, ySep, zSep]))
    else:
        axisVector = np.array([xSep, ySep, zSep])
    return axisVector


def normaliseVec(vector):
    '''This function normalises an input vector to unit magnitude'''
    return vector/float(np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2))


def findIndex(string, character):
    '''This function returns the locations of an inputted character in an inputted string'''
    index = 0
    locations = []
    while index < len(string):
        if string[index] == character:
            locations.append(index)
        index += 1
    if len(locations) == 0:
        return None
    return locations


def getRotationMatrix(vector1, vector2):
    '''This function returns the rotation matrix around the origin that maps vector1 to vector 2'''
    crossProduct = np.cross(vector1, vector2)
    sinAngle = np.sqrt(((crossProduct[0]**2) + ((crossProduct[1])**2) + ((crossProduct[2])**2)))
    cosAngle = np.dot(vector1, vector2)
    skewMatrix = np.matrix([[0, -crossProduct[2], crossProduct[1]], [crossProduct[2], 0, -crossProduct[0]], [-crossProduct[1], crossProduct[0], 0]])
    skewMatrixSquared = skewMatrix * skewMatrix
    rotMatrix = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) + skewMatrix + skewMatrixSquared*((1 - cosAngle)/(sinAngle**2))
    return rotMatrix


def parallelSort(list1, list2):
    '''This function sorts a pair of lists by the first list in ascending order (for example, atom mass and corresponding position can be input, sorted by ascending mass, and the two lists output, where the mass[atom_i] still corresponds to position[atom_i]'''
    data = zip(list1, list2)
    data.sort()
    list1, list2 = map(lambda t: list(t), zip(*data))
    return list1, list2


def writeCSV(dataX, dataY, name):
    '''Appends a CSV file with X and Y Data'''
    filename = './'+name+'.csv'
    document = csv.writer(open(filename, 'a+'), delimiter = ',')
    document.writerow([dataX, dataY])


def rotationMatrix(vector1, vector2):
    '''A function to return the rotation matrix around the origin that maps vector1 to vector 2'''
    crossProduct = np.cross(vector1, vector2)
    sinAngle = np.sqrt(((crossProduct[0]**2) + ((crossProduct[1])**2) + ((crossProduct[2])**2)))
    cosAngle = np.dot(vector1, vector2)
    skewMatrix = np.matrix([[0, -crossProduct[2], crossProduct[1]], [crossProduct[2], 0, -crossProduct[0]], [-crossProduct[1], crossProduct[0], 0]])
    skewMatrixSquared = skewMatrix * skewMatrix
    rotMatrix = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) + skewMatrix + skewMatrixSquared*((1 - cosAngle)/(sinAngle**2))
    return rotMatrix


def addUnwrappedPositions(inputDictionary):
    '''This function takes a runHoomd.py input dictionary and updates the 'unwrapped_position' key based on the values of the 'position' and 'image' keys'''
    simulationDimensions = [inputDictionary['lx'], inputDictionary['ly'], inputDictionary['lz']]
    inputDictionary['unwrapped_position'] = [0]*len(inputDictionary['position'])
    for i in range(len(inputDictionary['position'])):
        position = inputDictionary['position'][i]
        if len(inputDictionary['image']) > 0:
            image = inputDictionary['image'][i]
        else:
            image = [0, 0, 0]
        unwrappedPosition = []
        for axis in range(len(image)):
            unwrappedPosition.append((image[axis]*simulationDimensions[axis])+position[axis])
        # print "Original, Wrapped Position =", position, image
        # print "New, Unwrapped Position =", unwrappedPosition
        # raw_input("Press Return to continue...")
        inputDictionary['unwrapped_position'][i] = unwrappedPosition
    return inputDictionary


# def addWrappedPositions(inputDictionary):
#     inputDictionary['position'] = inputDictionary['unwrapped_position']
#     return inputDictionary


def addWrappedPositions(inputDictionary):
    '''This function takes a runHoomd.py input dictionary and updates the 'position' and 'image' keys based on the values of the 'unwrapped_position' key'''
    simulationDimensions = [inputDictionary['lx'], inputDictionary['ly'], inputDictionary['lz']]
    inputDictionary['position'] = [0]*len(inputDictionary['unwrapped_position'])
    inputDictionary['image'] = [0]*len(inputDictionary['unwrapped_position'])
    for atomID in range(len(inputDictionary['unwrapped_position'])):
        position = copy.deepcopy(inputDictionary['unwrapped_position'][atomID])
        imageCoords = [0, 0, 0]
        for axis in range(len(position)):
            if position[axis] > (simulationDimensions[axis]/2.0):
                while position[axis] > (simulationDimensions[axis]/2.0):
                    imageCoords[axis] += 1
                    position[axis] -= simulationDimensions[axis]
            elif position[axis] < -(simulationDimensions[axis]/2.0):
                while position[axis] < -(simulationDimensions[axis]/2.0):
                    imageCoords[axis] -= 1
                    position[axis] += simulationDimensions[axis]
        inputDictionary['position'][atomID] = position
        inputDictionary['image'][atomID] = imageCoords
    return inputDictionary


def addMasses(inputDictionary):
    '''This function takes a runHoomd.py input dictionary and updates the 'mass' key based on the values of the 'type' key. Note that more hardcoding is required to add aditional atom types'''
    inputDictionary['mass'] = [1.0]*len(inputDictionary['type'])
    for atomID in range(len(inputDictionary['type'])):
        if 'H' in inputDictionary['type'][atomID]:
            inputDictionary['mass'][atomID] = 1.00794
        elif 'C' in inputDictionary['type'][atomID]:
            inputDictionary['mass'][atomID] = 12.0107
        elif 'S' in inputDictionary['type'][atomID]:
            inputDictionary['mass'][atomID] = 32.0660
    return inputDictionary


def addDiameters(inputDictionary):
    '''This function takes a runHoomd.py input dictionary and updates the 'diameter' key based on the values of the 'type' key. Note that more hardcoding is required to add aditional atom types'''
    inputDictionary['diameter'] = [1.0]*len(inputDictionary['type'])
    for atomID in range(len(inputDictionary['type'])):
        if 'H' in inputDictionary['type'][atomID]:
            inputDictionary['diameter'][atomID] = 0.53
        elif 'C' in inputDictionary['type'][atomID]:
            inputDictionary['diameter'][atomID] = 0.67
        elif 'S' in inputDictionary['type'][atomID]:
            inputDictionary['diameter'][atomID] = 0.88
    return inputDictionary


def addTerminatingHydrogens(inputDictionary):
    '''This function takes a runHoomd.py input dictionary, determines the terminating monomers in the system, and creates new hydrogen atoms at a sensible position, bonded to the terminating carbon atoms of the chain'''
    # Find terminating carbon on the molecule
    connectingCarbons = {}
    terminatingCarbonIDs = []
    terminatingHydrogenIDs = []
    for bond in inputDictionary['bond']:
        if 'C1' in bond[0]:
            # Should pick up C1 and C10
            if 'C1' in inputDictionary['type'][bond[1]]:
                if bond[1] not in connectingCarbons:
                    connectingCarbons[bond[1]] = [1]
                else:
                    connectingCarbons[bond[1]][0] += 1
                if 'S1' in inputDictionary['type'][bond[2]]:
                    connectingCarbons[bond[1]].append(bond[2])
            if 'C1' in inputDictionary['type'][bond[2]]:
                if bond[2] not in connectingCarbons:
                    connectingCarbons[bond[2]] = [1]
                else:
                    connectingCarbons[bond[2]][0] += 1
                if 'S1' in inputDictionary['type'][bond[1]]:
                    connectingCarbons[bond[2]].append(bond[1])
    for connectingCarbonID in connectingCarbons:
        # Terminating Hydrogens are needed when there are only 2 bonds
        # (Hardcoded for P3HT)
        if connectingCarbons[connectingCarbonID][0] == 2:
            terminatingCarbonIDs.append([connectingCarbonID, connectingCarbons[connectingCarbonID][1]])
    # Now add the hydrogens
    for terminatingCarbon in terminatingCarbonIDs:
        # Get Hydrogen Position
        # Put the Hydrogen 1 unit (angstroems currently) away from the carbon, in the opposite direction to the sulfur.
        carbonPosition = inputDictionary['position'][terminatingCarbon[0]]
        sulfurPosition = inputDictionary['position'][terminatingCarbon[1]]
        sulfurAxis = findAxis(sulfurPosition, carbonPosition)
        hydrogenPosition = list(np.array(carbonPosition) + sulfurAxis)
        # Update the dictionary with all of the associated information for the hydrogen
        inputDictionary['position'].append(hydrogenPosition)
        inputDictionary['image'].append(inputDictionary['image'][terminatingCarbon[0]])
        # Add the unwrapped_position
        unwrappedPosition = []
        image = inputDictionary['image'][terminatingCarbon[0]]
        simulationDimensions = [inputDictionary['lx'], inputDictionary['ly'], inputDictionary['lz']]
        for axis in range(len(image)):
            unwrappedPosition.append((image[axis]*simulationDimensions[axis])+hydrogenPosition[axis])
        inputDictionary['unwrapped_position'].append(unwrappedPosition)
        inputDictionary['mass'].append(1.00794)
        inputDictionary['diameter'].append(0.53)
        inputDictionary['type'].append('H1')
        inputDictionary['velocity'].append(inputDictionary['velocity'][terminatingCarbon[0]])
        inputDictionary['body'].append(inputDictionary['body'][terminatingCarbon[0]])
        inputDictionary['charge'].append(inputDictionary['charge'][terminatingCarbon[0]])
        # Add the bond information (angles and dihedrals are unnecessary)
        newBond = [inputDictionary['type'][terminatingCarbon[0]]+'-H1', terminatingCarbon[0], len(inputDictionary['type'])-1]
        inputDictionary['bond'].append(newBond)
        # Finally, update the number of atoms in the system (we just added one!)
        terminatingHydrogenIDs.append(inputDictionary['natoms'])
        inputDictionary['natoms'] += 1
    return inputDictionary, terminatingHydrogenIDs


def loadMorphologyXML(xmlPath, sigma=1.0):
    # XML has SimDims as <box
    # Positions as <position and <image
    # Velocities as <velocity
    # Mass as <mass
    # Diameters as <diameter
    # Type as <type
    # "Body" as <body
    # Bonds as <bond, with each line as bondA, bondB, etc.
    # Angles as <angle, with each angle as angleA, angleB, etc.
    # Dihedral as <dihedral
    # Improper as <improper (usually none in xml)
    # Charge as <charge
    AtomDictionary = {'position':[], 'image':[], 'velocity':[], 'mass':[], 'diameter':[], 'type':[], 'body':[], 'bond':[], 'angle':[], 'dihedral':[], 'improper':[], 'charge':[]}
    record = False
    with open(xmlPath, 'r') as xmlFile:
        xmlData = xmlFile.readlines()
        for line in xmlData:
            if ('</' in line):
                record = False
            elif ('<configuration' in line) or ('<box' in line):
                # Get configuration data from this line (timestep, natoms etc)
                splitLine = line.split(' ')
                for i in range(1,len(splitLine)):
                    equalsLoc = findIndex(splitLine[i], '=')
                    if equalsLoc == None:
                        # Skip any elements without equals
                        continue
                    quotationLoc = findIndex(splitLine[i], '"')
                    if ('.' in splitLine[i][quotationLoc[0]+1:quotationLoc[1]]):
                        # Catch float in the value (excludes the = and quotation marks)
                        if ('<box' in line):
                            AtomDictionary[splitLine[i][:equalsLoc[0]]] = float(splitLine[i][quotationLoc[0]+1:quotationLoc[1]])*sigma
                        else:
                            AtomDictionary[splitLine[i][:equalsLoc[0]]] = float(splitLine[i][quotationLoc[0]+1:quotationLoc[1]])
                    else:
                        if ('<box' in line):
                            AtomDictionary[splitLine[i][:equalsLoc[0]]] = int(splitLine[i][quotationLoc[0]+1:quotationLoc[1]])*sigma
                        else:
                            AtomDictionary[splitLine[i][:equalsLoc[0]]] = int(splitLine[i][quotationLoc[0]+1:quotationLoc[1]])
            elif ('<position' in line):
                record = True
                recordType = 'position'
                continue
            elif ('<image' in line):
                record = True
                recordType = 'image'
                continue
            elif ('<velocity' in line):
                record = True
                recordType = 'velocity'
                continue
            elif ('<mass' in line):
                record = True
                recordType = 'mass'
                continue
            elif ('<diameter' in line):
                record = True
                recordType = 'diameter'
                continue
            elif ('<type' in line):
                record = True
                recordType = 'type'
                continue
            elif ('<body' in line):
                record = True
                recordType = 'body'
                continue
            elif ('<bond' in line):
                record = True
                recordType = 'bond'
                continue
            elif ('<angle' in line):
                record = True
                recordType = 'angle'
                continue
            elif ('<dihedral' in line):
                record = True
                recordType = 'dihedral'
                continue
            elif ('<improper' in line):
                record = True
                recordType = 'improper'
                continue
            elif ('<charge' in line):
                record = True
                recordType = 'charge'
                continue
            # Now we know what the variable is, append it to the dictionary data
            if (record == True):
                if (recordType == 'position') or (recordType == 'velocity'):
                    # NOTE: VELOCITIES ARE NOT NORMALISED IN THE MORPHOLOGY FILE...DO THEY NEED TO BE SCALED BY SIGMA OR NOT? CURRENTLY THEY ARE.
                    # Write to dictionary as floats scaled by sigma
                    splitLine = line.split(' ')
                    # Remove the "\n"
                    splitLine[-1] = splitLine[-1][:-1]
                    if (len(splitLine) == 1):
                        AtomDictionary[recordType].append(float(splitLine[0])*sigma)
                        continue
                    for i in range(len(splitLine)):
                        splitLine[i] = float(splitLine[i])*sigma
                    AtomDictionary[recordType].append(splitLine)
                elif (recordType == 'mass') or (recordType == 'diameter') or (recordType == 'charge'):
                    # Write to dictionary as floats
                    splitLine = line.split(' ')
                    # Remove the "\n"
                    splitLine[-1] = splitLine[-1][:-1]
                    if (len(splitLine) == 1):
                        AtomDictionary[recordType].append(float(splitLine[0]))
                        continue
                    for i in range(len(splitLine)):
                        splitLine[i] = float(splitLine[i])
                    AtomDictionary[recordType].append(splitLine)
                elif (recordType == 'image') or (recordType == 'body'):
                    # Write to dictionary as int
                    splitLine = line.split(' ')
                    # Remove the "\n"
                    splitLine[-1] = splitLine[-1][:-1]
                    if (len(splitLine) == 1):
                        AtomDictionary[recordType].append(int(splitLine[0]))
                        continue
                    for i in range(len(splitLine)):
                        splitLine[i] = int(splitLine[i])
                    AtomDictionary[recordType].append(splitLine)
                elif (recordType == 'type'):
                    # Write to dictionary as str
                    splitLine = line.split(' ')
                    # Remove the "\n"
                    splitLine[-1] = splitLine[-1][:-1]
                    AtomDictionary[recordType].append(str(splitLine[0]))
                else:
                    #  (recordType == 'bond') or (recordType == 'angle') or (recordType == 'dihedral') or (recordType == 'improper')
                    # Write to dictionary as combination
                    splitLine = line.split(' ')
                    # Remove the "\n"
                    splitLine[-1] = splitLine[-1][:-1]
                    splitLine[0] = str(splitLine[0])
                    for i in range(1,len(splitLine)):
                        splitLine[i] = int(splitLine[i])
                    AtomDictionary[recordType].append(splitLine)
    return AtomDictionary


def writeMorphologyXML(inputDictionary, outputFile):
    # First, need to check the positions of the atoms to ensure that everything is correctly contained inside the box
    print "Checking wrapped positions before writing XML..."
    inputDictionary = checkWrappedPositions(inputDictionary)
    # inputDictionary['position'], inputDictionary['image'] = pbc.shift_pbc(inputDictionary['position'], [inputDictionary['lx'], inputDictionary['ly'], inputDictionary['lz']])
    # print inputDictionary['image'][:20]

    # raw_input('HALT')

    
    # Add Boiler Plate first
    linesToWrite = ['<?xml version="1.0" encoding="UTF-8"?>\n', '<hoomd_xml version="1.4">\n', '<configuration time_step="0" dimensions="3" natoms="'+str(inputDictionary['natoms'])+'" >\n', '<box lx="'+str(inputDictionary['lx'])+'" ly="'+str(inputDictionary['ly'])+'" lz="'+str(inputDictionary['lz'])+'" />\n']
    # Position
    linesToWrite.append('<position num="'+str(inputDictionary['natoms'])+'">\n')
    for positionData in inputDictionary['position']:
        linesToWrite.append(" ".join(str(coord) for coord in positionData)+'\n')
    linesToWrite.append('</position>\n')
    # Image
    linesToWrite.append('<image num="'+str(inputDictionary['natoms'])+'">\n')
    for imageData in inputDictionary['image']:
        linesToWrite.append(" ".join(str(coord) for coord in imageData)+'\n')
    linesToWrite.append('</image>\n')
    # Velocity
    linesToWrite.append('<velocity num="'+str(inputDictionary['natoms'])+'">\n')
    for velocityData in inputDictionary['velocity']:
        linesToWrite.append(" ".join(str(coord) for coord in velocityData)+'\n')
    linesToWrite.append('</velocity>\n')
    # Mass
    linesToWrite.append('<mass num="'+str(inputDictionary['natoms'])+'">\n')
    for massData in inputDictionary['mass']:
        linesToWrite.append(str(massData)+'\n')
    linesToWrite.append('</mass>\n')
    # Diameter
    linesToWrite.append('<diameter num="'+str(inputDictionary['natoms'])+'">\n')
    for diameterData in inputDictionary['diameter']:
        linesToWrite.append(str(diameterData)+'\n')
    linesToWrite.append('</diameter>\n')
    # Type
    linesToWrite.append('<type num="'+str(inputDictionary['natoms'])+'">\n')
    for typeData in inputDictionary['type']:
        linesToWrite.append(str(typeData)+'\n')
    linesToWrite.append('</type>\n')
    # Body
    linesToWrite.append('<body num="'+str(inputDictionary['natoms'])+'">\n')
    for bodyData in inputDictionary['body']:
        linesToWrite.append(str(bodyData)+'\n')
    linesToWrite.append('</body>\n')
    # Bond
    linesToWrite.append('<bond num="'+str(len(inputDictionary['bond']))+'">\n')
    for bondData in inputDictionary['bond']:
        linesToWrite.append(" ".join(str(coord) for coord in bondData)+'\n')
    linesToWrite.append('</bond>\n')
    # Angle
    linesToWrite.append('<angle num="'+str(len(inputDictionary['angle']))+'">\n')
    for angleData in inputDictionary['angle']:
        linesToWrite.append(" ".join(str(coord) for coord in angleData)+'\n')
    linesToWrite.append('</angle>\n')
    # Dihedral
    linesToWrite.append('<dihedral num="'+str(len(inputDictionary['dihedral']))+'">\n')
    for dihedralData in inputDictionary['dihedral']:
        linesToWrite.append(" ".join(str(coord) for coord in dihedralData)+'\n')
    linesToWrite.append('</dihedral>\n')
    # Improper
    linesToWrite.append('<improper num="'+str(len(inputDictionary['improper']))+'">\n')
    for improperData in inputDictionary['improper']:
        linesToWrite.append(" ".join(str(coord) for coord in improperData)+'\n')
    linesToWrite.append('</improper>\n')
    # Charge
    linesToWrite.append('<charge num="'+str(inputDictionary['natoms'])+'">\n')
    for chargeData in inputDictionary['charge']:
        linesToWrite.append(str(chargeData)+'\n')
    linesToWrite.append('</charge>\n')
    linesToWrite.append('</configuration>\n')
    linesToWrite.append('</hoomd_xml>\n')
    with open(outputFile, 'w+') as xmlFile:
        xmlFile.writelines(linesToWrite)


def writePOSCARFile(inputDict, outputFile):
    '''This function takes an input dictionary and converts it to a POSCAR for use in DFT calculations'''
    # This POSCAR is ordered as C, S, H for Izaak.
    linesToWrite = []
    atomsByType = [[],[],[]] # C, S, H
    for atomID in range(len(inputDict['type'])):
        atomType = inputDict['type'][atomID][0]
        if atomType == 'C':
            atomsByType[0].append(atomID)
        elif atomType == 'S':
            atomsByType[1].append(atomID)
        elif atomType == 'H':
            atomsByType[2].append(atomID)
    
    # linesToWrite = []
    # typeList = []
    # freqList = []
    # previousType = None
    # numberOfTypes = 0
    # for atomID in range(len(inputDict['type'])):
    #     atomType = inputDict['type'][atomID][0]
    #     if atomType != previousType:
    #         if previousType != None:
    #             typeList.append(previousType)
    #             freqList.append(numberOfTypes)
    #         previousType = atomType
    #         numberOfTypes = 1
    #     else:
    #         numberOfTypes += 1
    # # Now have to add the final lot of atoms:
    # typeList.append(previousType)
    # freqList.append(numberOfTypes)
    # Line 1 = CommentLine
    slashLocs = findIndex(outputFile, '/')
    linesToWrite.append(str(outputFile[slashLocs[-3]+1:slashLocs[-2]+1])+str(outputFile[slashLocs[-1]+1:]).replace('.POSCAR', '')+' VASP input file.\n')
    # Line 2 = Scale Factor
    linesToWrite.append('1.000000000000\n')
    # Lines 3-5 = Box Dimensions
    boxDims = []
    for key in ['lx', 'ly', 'lz']:
        boxDims.append(inputDict[key])
    boxDimsMatrix = np.diag(np.array(boxDims))
    for row in boxDimsMatrix:
        boxRow = ''
        for element in row:
            boxRow += "{:22.15f}".format(element)
        linesToWrite.append(boxRow+'\n')
    # Line 6 = Atom Types
    # linesToWrite.append(' '.join(typeList)+'\n')
    linesToWrite.append('C S H \n')
    # Line 7 = Frequency of Types
    # linesToWrite.append(' '.join(map(str, freqList))+'\n')
    linesToWrite.append(str(len(atomsByType[0]))+' '+str(len(atomsByType[1]))+' '+str(len(atomsByType[2]))+'\n')
    # Line 8 = 'Cartesian'
    linesToWrite.append('Cartesian\n')
    # Lines 9+ = Positions
    # Positions are not set to be origin in the middle, origin is bottom left corner. As such, we need to add L/2 to each coordinate
    writeOrder = []
    for atomType in atomsByType:
        writeOrder += atomType
    # for position in inputDict['position']:
    #     coordinates = ''
    #     for axis in range(3):
    #         coordinates += "{:22.15f}".format(position[axis]+(boxDims[axis]/2.))
    #     linesToWrite.append(coordinates+'\n')
    for atomID in writeOrder:
        coordinates = ''
        for axis in range(3):
            coordinates += "{:22.15f}".format(inputDict['position'][atomID][axis]+(boxDims[axis]/2.))
        linesToWrite.append(coordinates+'\n')
    with open(outputFile, 'w+') as POSCARFile:
        POSCARFile.writelines(linesToWrite)
    with open(outputFile.replace('POSCAR', 'pickle'), 'w+') as bondPickle:
        pickle.dump(inputDict['bond'], bondPickle)
    print "POSCAR data written to", str(outputFile)+". Bond data written to", str(outputFile.replace('POSCAR', 'pickle'))+"."

        
def writeXYZFile(inputDict, outputFile):
    '''This function takes an input dictionary and converts it to an XYZ for use in DFT calculations'''
    # First line is atom numbers, second line is boiler plate
    rowsToWrite = [str(inputDict['natoms'])+'\n', 'XYZ file generated from XML using helperFunctions.XMLToXYZ\n']
    # Format of xyz is Type, X Pos, Y Pos, Z Pos
    for atomID in range(len(inputDict['type'])):
        # Note, this will break for atoms that have two-letter symbols (e.g. Al, Ca etc.)
        atomType = inputDict['type'][atomID][0]
        while len(atomType) < 10:
            atomType += ' '
        atomX = str(inputDict['position'][atomID][0])
        while len(atomX) < 20:
            atomX += ' '
        atomY = str(inputDict['position'][atomID][1])
        while len(atomY) < 20:
            atomY += ' '
        atomZ = str(inputDict['position'][atomID][2])
        lineToWrite = atomType+atomX+atomY+atomZ+'\n'
        rowsToWrite.append(lineToWrite)
    with open(outputFile, 'w+') as xyzFile:
        xyzFile.writelines(rowsToWrite)
    print "XYZ data written to", str(outputFile)+"."


def createSlurmSubmissionScript(outputDir, runName, mode):
    '''This function creates a slurm submission script for Kestrel from a template file sample.sh'''
    queue = 'batch'
    jobName = str(runName)
    outputFile = outputDir+'/'+str(runName)[:-4]+'.o'
    with open(os.getcwd()+'/templates/sample.sh', 'r') as template:
        templateLines = template.readlines()
    for lineNo in range(len(templateLines)):
        if '-p batch' in templateLines[lineNo]:
            # This is queue select
            templateLines[lineNo] = templateLines[lineNo].replace('batch', queue)
        elif '-J JOBNAME' in templateLines[lineNo]:
            # This is job name
            templateLines[lineNo] = templateLines[lineNo].replace('JOBNAME', jobName)
        elif '-o OUTFILE' in templateLines[lineNo]:
            # This is outfile
            templateLines[lineNo] = templateLines[lineNo].replace('OUTFILE', outputFile)
        elif '--mail-user' in templateLines[lineNo]:
            # E-mail address
            templateLines[lineNo] = templateLines[lineNo].replace('CHANGEME', 'mattyjones')
        elif '-t 12:00:00' in templateLines[lineNo]:
            # Wallclock time
            templateLines[lineNo] = templateLines[lineNo].replace('12:00:00', '01:00:00')
        elif 'cd /scratch/${USER}' in templateLines[lineNo]:
            templateLines[lineNo] = templateLines[lineNo].replace('/scratch/${USER}', os.getcwd())
        elif 'myfile.py' in templateLines[lineNo]:
            # This is actual execute line
            templateLines[lineNo] = templateLines[lineNo].replace('myfile.py', os.getcwd()+'/code/'+'runHoomd.py '+outputDir+'/'+runName)
            if mode == 'cpu':
                templateLines[lineNo] = templateLines[lineNo].replace('--mode=gpu --gpu=0', '--mode=cpu')
        # Finally, sort out the /scratch/ space and move the output files somewhere useful
    submissionScriptName = outputDir+'/'+runName[:-4]+'.sh'
    with open(submissionScriptName, 'w+') as submissionScript:
        submissionScript.writelines(templateLines)
    return submissionScriptName


def incrementAtomIDs(inputDictionary, increment):
    for bond in inputDictionary['bond']:
        bond[1] += increment
        bond[2] += increment
    for angle in inputDictionary['angle']:
        angle[1] += increment
        angle[2] += increment
        angle[3] += increment
    for dihedral in inputDictionary['dihedral']:
        dihedral[1] += increment
        dihedral[2] += increment
        dihedral[3] += increment
        dihedral[4] += increment
    for improper in inputDictionary['improper']:
        improper[1] += increment
        improper[2] += increment
        improper[3] += increment
        improper[4] += increment
    return inputDictionary


def scale(inputDictionary, scaleFactor):
    for ID, position in enumerate(inputDictionary['position']):
        #if ID == 24104:
            # print "Initial Position =", inputDictionary['position'][ID], inputDictionary['image'][ID]
        inputDictionary['position'][ID] = list(scaleFactor*np.array(position))
        #if ID == 24104:
            # print "Scaled Position =", inputDictionary['position'][ID], inputDictionary['image'][ID]
    for element in ['lx', 'ly', 'lz']:
        if element in inputDictionary:
            inputDictionary[element] *= scaleFactor
	ipos = np.array(inputDictionary['position'])
    return inputDictionary


def centre(inputDictionary, centreOfMass):
    COM = np.array(centreOfMass)
    numberOfAtomsMoved = 0
    for index, position in enumerate(inputDictionary['position']):
        distanceMoved0 = calculateSeparation(position, COM)
        inputDictionary['position'][index] = list(position-COM)
        distanceMoved = calculateSeparation(position, COM)
        # print "Distance Moved by atom", index, "=", distanceMoved0, distanceMoved
        numberOfAtomsMoved += 1
    return inputDictionary


def checkWrappedPositions(inputDictionary):
    atomPositions = np.array(inputDictionary['position'])
    atomImages = np.array(inputDictionary['image'])
    xhi = inputDictionary['lx']/2.0
    xlo = -inputDictionary['lx']/2.0
    yhi = inputDictionary['ly']/2.0
    ylo = -inputDictionary['ly']/2.0
    zhi = inputDictionary['lz']/2.0
    zlo = -inputDictionary['lz']/2.0
    # tp=pbc.plain_pbc(atomPositions,(inputDictionary['lx'],inputDictionary['ly'],inputDictionary['lz']) )
    # tp=pbc.plain_pbc(tp,(inputDictionary['lx'],inputDictionary['ly'],inputDictionary['lz']) )
    # tp=pbc.plain_pbc(tp,(inputDictionary['lx'],inputDictionary['ly'],inputDictionary['lz']) )
    # tp=pbc.plain_pbc(tp,(inputDictionary['lx'],inputDictionary['ly'],inputDictionary['lz']) )
    # tp=pbc.plain_pbc(tp,(inputDictionary['lx'],inputDictionary['ly'],inputDictionary['lz']) )
    # tp=pbc.plain_pbc(tp,(inputDictionary['lx'],inputDictionary['ly'],inputDictionary['lz']) )
    # tp=pbc.plain_pbc(tp,(inputDictionary['lx'],inputDictionary['ly'],inputDictionary['lz']) )
    # tp=pbc.plain_pbc(tp,(inputDictionary['lx'],inputDictionary['ly'],inputDictionary['lz']) )
    
    for atomID in range(len(atomPositions)):
        while atomPositions[atomID][0] > xhi:
            atomPositions[atomID][0] -= inputDictionary['lx']
            atomImages[atomID][0] += 1
        while atomPositions[atomID][0] < xlo:
            atomPositions[atomID][0] += inputDictionary['lx']
            atomImages[atomID][0] -= 1
        while atomPositions[atomID][1] > yhi:
            atomPositions[atomID][1] -= inputDictionary['ly']
            atomImages[atomID][1] += 1
        while atomPositions[atomID][1] < ylo:
            atomPositions[atomID][1] += inputDictionary['ly']
            atomImages[atomID][1] -= 1
        while atomPositions[atomID][2] > zhi:
            atomPositions[atomID][2] -= inputDictionary['lz']
            atomImages[atomID][2] += 1
        while atomPositions[atomID][2] < zlo:
            atomPositions[atomID][2] += inputDictionary['lz']
            atomImages[atomID][2] -= 1
    inputDictionary['position'] = list(atomPositions)
    inputDictionary['image'] = list(atomImages)
    # print np.sum(np.absolute(atomPositions-tp) > 0.)
    return inputDictionary


def alignMolecule(inputDictionary, vectorToAlignTo):
    '''This function rotates a molecule such that the vector between the first and last sulfur atoms in the chain (assumed
    to be the backbone vector) is mapped to vectorToAlignTo'''
    sulfurAtomIDs = []
    for atomIndex, atomType in enumerate(inputDictionary['type']):
        if atomType[0] == 'S':
            sulfurAtomIDs.append(atomIndex)
    sulfurAtomIDs.sort()
    chainOrientationVector = findAxis(inputDictionary['position'][sulfurAtomIDs[0]], inputDictionary['position'][sulfurAtomIDs[-1]])
    vectorToAlignTo = np.array(vectorToAlignTo)
    rotationMatrix = getRotationMatrix(chainOrientationVector, vectorToAlignTo)
    for atomID, pos in enumerate(inputDictionary['position']):
        positionArray = np.copy(pos)
        rotatedPosition = np.transpose(rotationMatrix*np.transpose(np.matrix(positionArray)))
        inputDictionary['position'][atomID] = [rotatedPosition[0,0], rotatedPosition[0,1], rotatedPosition[0,2]]
    return inputDictionary


def cellSearchBonds(moleculeDict):
    '''This function finds the bonds in the system based on the proximity of atoms to their neighbours'''
    raise SystemError("THIS FUNCTION DOES NOT WORK AND IT'S NONE-TRIVIAL TO IMPLEMENT")
    moleculeDict['neighbourCell'] = []
    maximumBondLength = 1.6# Bond length in angstroems
    atomIDs = np.arange(len(moleculeDict['position']))
    for coordinates in moleculeDict['position']:
        cellLocation = np.copy(coordinates)
        moleculeDict['neighbourCell'].append(map(int, np.round(cellLocation/maximumBondLength)))
        print coordinates, moleculeDict['neighbourCell'][-1]
    neighbourCells = np.copy(moleculeDict['neighbourCell'])
    
    print calculateSeparation(moleculeDict['position'][0], moleculeDict['position'][1]), calculateSeparation(moleculeDict['position'][1], moleculeDict['position'][2])
    parallelSort(neighbourCells, atomIDs)
    for i in range(len(atomIDs)):
        print atomIDs[i], neighbourCells[i]
        

def getAAIDsByMolecule(CGtoAAIDs):
    '''This function extracts the molecule AAIDs given a dictionary CGtoAAIDs which describes the mapping of all atom particles to each CG site'''
    moleculeAAIDs = []
    for moleculeID, CGtoAAIDDict in enumerate(CGtoAAIDs):
        moleculeAAIDs.append([])
        for dictionaryValue in CGtoAAIDs[moleculeID].values():
            moleculeAAIDs[-1] += dictionaryValue[1]
        moleculeAAIDs[-1].sort()
    return moleculeAAIDs


def getsScale(outputDir, morphologyName):
    morphologyFiles = os.listdir(outputDir+'/morphology')
    for fileName in morphologyFiles:
        if 'scaled' in fileName:
            scaledXMLName = fileName
            break
    underscoreLocs = findIndex(scaledXMLName, '_')
    inverseScaleFactor = scaledXMLName[underscoreLocs[-2]+1:underscoreLocs[-1]]
    return float(inverseScaleFactor)


def loadDict(masterDict, moleculeIDs, bondPickleName):
    '''This function generates a molecule dictionary by picking the relevant data from a masterDict using a list of atomIDs given by moleculeIDs'''
    moleculeDict = {'position':[], 'unwrapped_position':[], 'type':[], 'diameter':[], 'image':[], 'charge':[], 'mass':[], 'velocity':[]}
    # First get atom-specific properties
    for atomID in moleculeIDs:
        for key in moleculeDict.keys():
            moleculeDict[key].append(masterDict[key][atomID])
    # Then add in the simulation properties
    for key in ['lx', 'ly', 'lz', 'xy', 'xz', 'yz', 'dimensions']:
        moleculeDict[key] = masterDict[key]
    # Then load the relevant bonds
    with open(bondPickleName, 'r') as bondPickle:
        moleculeDict['bond'] = pickle.load(bondPickle)
    # Now need to unwrap the coordinates
    # moleculeDict = addUnwrappedPositions(moleculeDict)
    # # Set the unwrapped coordinates to the default 'position' (saves time on rewriting some analyseMolecules functions and shouldn't affect anything)
    # moleculeDict['position'] = moleculeDict['unwrapped_position']
    moleculeCOM = calcCOM(moleculeDict['position'], moleculeDict['mass'])
    return moleculeDict
    
        
def loadPoscar(inputFilePath):
    '''This function loads a poscar file located at inputFilePath, and creates a dictionary of the atomic types and positions.
    It also loads the pickle file containing the bond information and adds it to the dictionary.'''
    moleculeDict = {'position':[], 'type':[]}
    with open(inputFilePath, 'r') as poscarFile:
        poscarData = poscarFile.readlines()
    simDims = []
    for unitCellLine in poscarData[2:5]:
        simDims.append([])
        for coordinate in unitCellLine[:-1].split(' '):
            if len(coordinate) > 0:
                simDims[-1].append(float(coordinate))
    moleculeDict['lx'] = simDims[0][0]
    moleculeDict['ly'] = simDims[1][1]
    moleculeDict['lz'] = simDims[2][2]
    simBoxDims = ['lx', 'ly', 'lz']
    typeList = poscarData[5].split('\n')[0].split(' ')
    freqList = map(int, poscarData[6].split('\n')[0].split(' '))
    for i in range(len(typeList)):
        if len(typeList[i]) != 0: # Catch for the extra space I had to put in to make VMD behave properly
            moleculeDict['type'] += [typeList[i]]*freqList[i]
    for atomCoords in poscarData[8:]:
        moleculeDict['position'].append([])
        for coordinate in atomCoords.split('\n')[0].split(' '):
            coordinatesToWrite = []
            if len(coordinate) > 0:
                coordinatesToWrite.append(float(coordinate))
            for i in range(len(coordinatesToWrite)):
                moleculeDict['position'][-1].append(coordinatesToWrite[i]-(moleculeDict[simBoxDims[i]]/2.0))
    with open(inputFilePath.replace('POSCAR', 'pickle'), 'r') as bondPickle:
        moleculeDict['bond'] = pickle.load(bondPickle)
    moleculeDict = addMasses(moleculeDict)
    return moleculeDict


def checkORCAFileStructure(outputDir):
    '''This function checks that the correct directories are in place for the ORCA transfer-integral calculation'''
    morphologyDirList = os.listdir(outputDir)
    if 'chromophores' not in morphologyDirList:
        print "Making /chromophores directory..."
        os.makedirs(outputDir+'/chromophores')
        print "Making /inputORCA directory..."
        os.makedirs(outputDir+'/chromophores/inputORCA')
        os.makedirs(outputDir+'/chromophores/inputORCA/single')
        os.makedirs(outputDir+'/chromophores/inputORCA/pair')
        print "Making /outputORCA directory..."
        os.makedirs(outputDir+'/chromophores/outputORCA')
    else:
        chromophoresDirList = os.listdir(outputDir+'/chromophores')
        if 'inputORCA' not in chromophoresDirList:
            print "Making /inputORCA directory..."
            os.makedirs(outputDir+'/chromophores/inputORCA')
            os.makedirs(outputDir+'/chromophores/inputORCA/single')
            os.makedirs(outputDir+'/chromophores/inputORCA/pair')
        if 'outputORCA' not in chromophoresDirList:
            print "Making /outputORCA directory..."
            os.makedirs(outputDir+'/chromophores/outputORCA')


def writeORCAInp(inputDictList, outputDir, mode):
    '''This function loads the ORCA input template and creates the segment pair ORCA inputs for this morphology, for running later'''
    chromophore1 = inputDictList[0]
    chromo1Name = str(chromophore1['realChromoID'])
    while len(chromo1Name) < 4:
        chromo1Name = '0'+chromo1Name
    if mode == 'pair':
        chromophore2 = inputDictList[1]
        # First check that the file doesn't already exist
        chromo2Name = str(chromophore2['realChromoID'])
        while len(chromo2Name) < 4:
            chromo2Name = '0'+chromo2Name
        ORCAFileName = 'chromo'+chromo1Name+'_chromo'+chromo2Name+'.inp'
    elif mode == 'single':
        ORCAFileName = 'chromo'+chromo1Name+'.inp'
    # Check by opening the file - saves time on regenerating the os.listdirs list for many thousands of files
    try:
        with open(outputDir+'/chromophores/inputORCA/'+mode+'/'+ORCAFileName, 'r') as testFile:
            fileExists = True
    except IOError:
        fileExists = False
    inputFileName = outputDir+'/chromophores/inputORCA/'+mode+'/'+ORCAFileName
    if fileExists == True:
        print "\n"
        print "File", ORCAFileName, "already exists, skipping..."
        #return
        print "Creating file anyway to check that they are the same"
        inputFileName = inputFileName.replace('.inp', '_2.inp')
    if mode == 'pair':
        # Centre the dimer pair at the origin
        COM = calcCOM(chromophore1['position']+chromophore2['position'], chromophore1['mass']+chromophore2['mass'])
    elif mode == 'single':
        # Centre the chromophore at the origin
        COM = calcCOM(chromophore1['position'], chromophore1['mass'])
    chromophore1 = centre(chromophore1, COM)
    if mode == 'pair':
        chromophore2 = centre(chromophore2, COM)
    # Now write the file
    with open(os.getcwd()+'/templates/template.inp', 'r') as templateFile:
        inpFileLines = templateFile.readlines()
    linesToWrite = []
    for atomID, atomCoords in enumerate(chromophore1['position']):
        thisAtomData = ' '+chromophore1['type'][atomID][0]+' '+' '.join(map(str, atomCoords))+'\n'
        linesToWrite.append(thisAtomData)
    if mode == 'pair':
        for atomID, atomCoords in enumerate(chromophore2['position']):
            thisAtomData = ' '+chromophore2['type'][atomID][0]+' '+' '.join(map(str, atomCoords))+'\n'
            linesToWrite.append(thisAtomData)
    inpFileLines[-1:-1] = linesToWrite
    with open(inputFileName, 'w+') as ORCAInputFile:
        ORCAInputFile.writelines(inpFileLines)
    print "ORCA input file written to", inputFileName, "\r",
    if fileExists == True:
        raw_input("Hit return to continue...")

        
def getORCAJobs(inputDir):
    ORCAFileList = os.listdir(inputDir)
    ORCAFilesToRun = []
    for fileName in ORCAFileList:
        if '.inp' in fileName:
            ORCAFilesToRun.append(inputDir+'/'+fileName)
    ORCAFilesToRun.sort()
    try:
        procIDs = os.environ.get('SLURM_GTIDS').split(',')
    except AttributeError:
        # Was not loaded using SLURM, so use all physical processors
        procIDs = list(np.arange(mp.cpu_count()))
    jobsList = [ORCAFilesToRun[i:i+(len(ORCAFilesToRun)/len(procIDs))] for i in xrange(0, len(ORCAFilesToRun), int(np.ceil(len(ORCAFilesToRun)/float(len(procIDs)))))]
    return procIDs, jobsList
