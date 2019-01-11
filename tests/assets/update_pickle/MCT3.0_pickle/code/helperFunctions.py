import numpy as np
import copy
import os
import sys
import pickle
import multiprocessing as mp
import csv
import xml.etree.cElementTree as ET
import time as T


# UNIVERSAL CONSTANTS, DO NOT CHANGE!
elementaryCharge = 1.60217657E-19 # C
kB = 1.3806488E-23 # m^{2} kg s^{-2} K^{-1}
hbar = 1.05457173E-34 # m^{2} kg s^{-1}


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
    atom1 = np.array(atom1)
    atom2 = np.array(atom2)
    return np.sqrt(np.sum((atom1 - atom2)**2))


def calcCOM(listOfPositions, listOfAtomTypes=None, listOfMasses=None):
    '''This function calculates the centre of mass of a collection of sites/atoms (listOfPositions) with corresponding type (listOfAtomTypes) or mass (listOfMasses)
    If listOfMasses is not specified, then listOfAtomTypes MUST be.'''
    massWeighted = np.array([0.0, 0.0, 0.0])
    if listOfMasses is None:
        listOfMasses = []
        for atomType in listOfAtomTypes:
            # Masses obtained from nist.gov, for the atoms we are likely to simulate the most.
            # Add in new atoms here if your molecule requires it!
            if ('BR' in atomType) or ('Br' in atomType) or ('br' in atomType):
                print("Br 79 being used as the preferred isotope, change in helperFunctions.calcCOM if not!")
                listOfMasses.append(78.918338)
            elif ('SI' in atomType) or ('Si' in atomType) or ('si' in atomType):
                listOfMasses.append(27.976926)
            elif ('C' in atomType) or ('c' in atomType):
                listOfMasses.append(12.000000)
            elif ('H' in atomType) or ('h' in atomType):
                listOfMasses.append(1.007825)
            elif ('S' in atomType) or ('s' in atomType):
                listOfMasses.append(31.972071)
            elif ('O' in atomType) or ('o' in atomType):
                listOfMasses.append(15.994914)
            elif ('N' in atomType) or ('n' in atomType):
                listOfMasses.append(14.003074)
            elif (atomType == 'D') or (atomType == 'A'):
                listOfMasses.append(1.0)
            else:
                raise SystemError("Unknown atomic mass", atomType, "please hardcode into helperFunctions.calcCOM.")
    totalMass = np.sum(listOfMasses)
    for atomID, position in enumerate(listOfPositions):
        for axis in range(3):
            massWeighted[axis] += position[axis] * listOfMasses[atomID]
    return massWeighted / float(totalMass)


def findAxis(atom1, atom2, normalise=True):
    '''This function determines the normalised vector from the location of atom1 to atom2. The positions can enter as lists or arrays, but are output as arrays'''
    xSep = atom2[0] - atom1[0]
    ySep = atom2[1] - atom1[1]
    zSep = atom2[2] - atom1[2]
    if normalise is True:
        axisVector = normaliseVec(np.array([xSep, ySep, zSep]))
    else:
        axisVector = np.array([xSep, ySep, zSep])
    return axisVector


def normaliseVec(vector):
    '''This function normalises an input vector to unit magnitude'''
    return vector / float(np.sqrt(np.sum(vector)**2))


def getRotationMatrix(vector1, vector2):
    '''This function returns the rotation matrix around the origin that maps vector1 to vector 2'''
    crossProduct = np.cross(vector1, vector2)
    sinAngle = np.sqrt(((crossProduct[0]**2) + ((crossProduct[1])**2) + ((crossProduct[2])**2)))
    cosAngle = np.dot(vector1, vector2)
    skewMatrix = np.matrix([[0, -crossProduct[2], crossProduct[1]], [crossProduct[2], 0, -crossProduct[0]], [-crossProduct[1], crossProduct[0], 0]])
    skewMatrixSquared = skewMatrix * skewMatrix
    rotMatrix = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) + skewMatrix + skewMatrixSquared * ((1 - cosAngle) / (sinAngle**2))
    return rotMatrix


def parallelSort(list1, list2):
    '''This function sorts a pair of lists by the first list in ascending order (for example, atom mass and corresponding position can be input, sorted by ascending mass, and the two lists output, where the mass[atom_i] still corresponds to position[atom_i]'''
    list1, list2 = zip(*sorted(zip(list1, list2)))
    return list1, list2


def writeCSV(fileName, data):
    '''Writes a CSV file given a 2D array `data' of arbitrary size'''
    with open(fileName, 'w+') as csvFile:
        document = csv.writer(csvFile, delimiter=',')
        for row in data:
            document.writerow(list(row))
    print("CSV written to", fileName)


def rotationMatrix(vector1, vector2):
    '''A function to return the rotation matrix around the origin that maps vector1 to vector 2'''
    crossProduct = np.cross(vector1, vector2)
    sinAngle = np.sqrt(((crossProduct[0]**2) + ((crossProduct[1])**2) + ((crossProduct[2])**2)))
    cosAngle = np.dot(vector1, vector2)
    skewMatrix = np.matrix([[0, -crossProduct[2], crossProduct[1]], [crossProduct[2], 0, -crossProduct[0]], [-crossProduct[1], crossProduct[0], 0]])
    skewMatrixSquared = skewMatrix * skewMatrix
    rotMatrix = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) + skewMatrix + skewMatrixSquared * ((1 - cosAngle) / (sinAngle**2))
    return rotMatrix


def addUnwrappedPositions(inputDictionary):
    '''This function takes a runHoomd.py input dictionary and updates the 'unwrapped_position' key based on the values of the 'position' and 'image' keys'''
    simulationDimensions = [inputDictionary['lx'], inputDictionary['ly'], inputDictionary['lz']]
    inputDictionary['unwrapped_position'] = [0] * len(inputDictionary['position'])
    for i in range(len(inputDictionary['position'])):
        position = inputDictionary['position'][i]
        if len(inputDictionary['image']) > 0:
            image = inputDictionary['image'][i]
        else:
            image = [0, 0, 0]
        unwrappedPosition = []
        for axis in range(len(image)):
            unwrappedPosition.append((image[axis] * simulationDimensions[axis]) + position[axis])
        inputDictionary['unwrapped_position'][i] = unwrappedPosition
    return inputDictionary


def replaceWrappedPositions(inputDictionary):
    '''This function takes a morphCT input dictionary and replaces the 'position' and 'image' keys with the 'unwrapped_position' key and '[0, 0, 0]' respectively.'''
    for atomID, unwrapped_position in enumerate(inputDictionary['unwrapped_position']):
        inputDictionary['position'][atomID] = unwrapped_position
        inputDictionary['image'][atomID] = [0, 0, 0]
    return inputDictionary


def addWrappedPositions(inputDictionary):
    '''This function takes a runHoomd.py input dictionary and updates the 'position' and 'image' keys based on the values of the 'unwrapped_position' key'''
    simulationDimensions = [inputDictionary['lx'], inputDictionary['ly'], inputDictionary['lz']]
    inputDictionary['position'] = [0] * len(inputDictionary['unwrapped_position'])
    inputDictionary['image'] = [0] * len(inputDictionary['unwrapped_position'])
    for atomID in range(len(inputDictionary['unwrapped_position'])):
        position = copy.deepcopy(inputDictionary['unwrapped_position'][atomID])
        imageCoords = [0, 0, 0]
        for axis in range(len(position)):
            if position[axis] > (simulationDimensions[axis] / 2.0):
                while position[axis] > (simulationDimensions[axis] / 2.0):
                    imageCoords[axis] += 1
                    position[axis] -= simulationDimensions[axis]
            elif position[axis] < -(simulationDimensions[axis] / 2.0):
                while position[axis] < -(simulationDimensions[axis] / 2.0):
                    imageCoords[axis] -= 1
                    position[axis] += simulationDimensions[axis]
        inputDictionary['position'][atomID] = position
        inputDictionary['image'][atomID] = imageCoords
    return inputDictionary


def addMasses(inputDictionary):
    '''This function takes a runHoomd.py input dictionary and updates the 'mass' key based on the values of the 'type' key. Note that more hardcoding is required to add aditional atom types'''
    inputDictionary['mass'] = [1.0] * len(inputDictionary['type'])
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
    inputDictionary['diameter'] = [1.0] * len(inputDictionary['type'])
    for atomID in range(len(inputDictionary['type'])):
        if 'H' in inputDictionary['type'][atomID]:
            inputDictionary['diameter'][atomID] = 0.53
        elif 'C' in inputDictionary['type'][atomID]:
            inputDictionary['diameter'][atomID] = 0.67
        elif 'S' in inputDictionary['type'][atomID]:
            inputDictionary['diameter'][atomID] = 0.88
    return inputDictionary


def getTerminatingPositions(currentAtomPosn, bondedAtomPositions, numberOfUnitsToAdd):
    # Given a currentAtomPosn and several bondedAtomPositions we can add numberOfUnitsToAdd different terminating units to the currentAtom through a series of geometric checks.
    # First get the vector to the average position of the bonded neighbours
    hydrogenPositions = []
    averagePositionOfBondedAtoms = np.array([0.0, 0.0, 0.0])
    for bondedAtomPosn in bondedAtomPositions:
        bondVector = np.array(bondedAtomPosn) - currentAtomPosn
        bondVector /= np.linalg.norm(bondVector)
        averagePositionOfBondedAtoms += bondVector
    [x, y, z] = currentAtomPosn + (-1.06 * (averagePositionOfBondedAtoms / np.linalg.norm(averagePositionOfBondedAtoms)))
    if numberOfUnitsToAdd == 1:
        # Easy, this is the perylene code
        # Simply reverse the bonded vector and make it the hydrogen position at a distance of 1.06 angstroems
        hydrogenPositions.append(np.array([x, y, z]))
    # Initial position for all hydrogens
    elif numberOfUnitsToAdd == 2:
        # As above (to get the right plane), but then rotated +(109.5/2) degrees and -(109.5/2) degrees around the bonding axis
        rotationAxis = np.array(bondedAtomPositions[0]) - np.array(bondedAtomPositions[-1])
        rotationAxis /= np.linalg.norm(rotationAxis)
        # Rotation matrix calculations from: http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
        # The array that describes the 3D rotation of (x, y, z) around the point (a, b, c) through
        # the unit axis <u, v, w> by the angle theta is given by:
        # [ (a(v^2 + w^2) - u(bv + cw - ux - vy - wz))(1 - cos(theta)) + x*cos(theta) + (-cv + bw - wy + vz)sin(theta),
        #   (b(u^2 + w^2) - v(au + cw - ux - vy - wz))(1 - cos(theta)) + y*cos(theta) + (cu - aw + wx - uz)sin(theta),
        #   (c(u^2 + v^2) - w(au + bv - ux - vy - wz))(1 - cos(theta)) + z*cos(theta) + (-bu + av - vx + uy)sin(theta) ]
        [a, b, c] = currentAtomPosn
        [u, v, w] = rotationAxis
        for theta in [(109.5 / 2.0) * (np.pi / 180.0), -(109.5 / 2.0) * (np.pi / 180.0)]:
            newPosition = np.array([(a * (v**2 + w**2) - u * ((b * v) + (c * w) - (u * x) - (v * y) - (w * z))) * (1 - np.cos(theta)) + (x * np.cos(theta)) + ((-(c * v) + (b * w) - (w * y) + (v * z)) * np.sin(theta)),
                           (b * (u**2 + w**2) - v * ((a * u) + (c * w) - (u * x) - (v * y) - (w * z))) * (1 - np.cos(theta)) + (y * np.cos(theta)) + (((c * u)  - (a * w) + (w * x) - (u * z)) * np.sin(theta)),
                           (c * (u**2 + v**2) - w * ((a * u) + (b * v) - (u * x) - (v * y) - (w * z))) * (1 - np.cos(theta)) + (z * np.cos(theta)) + ((-(b * u) + (a * v) - (v * x) + (u * y)) * np.sin(theta))])
            hydrogenPositions.append(newPosition)
    elif numberOfUnitsToAdd == 3:
        # As for one (to get the right side of the bonded atom), rotate the first one up by 70.5 (180 - 109.5) and then rotate around by 109.5 degrees for the other two
        # The first hydrogen can be rotated around any axis perpendicular to the only bond present
        axisToBond = currentAtomPosn - np.array(bondedAtomPositions[0])
        # Now find one of the set of vectors [i, j, k] perpendicular to this one so we can place the first hydrogen.
        # Do this by setting i = j = 1 and solve for k (given that currentAtomPosn[0]*i + currentAtomPosn[1]*j + currentAtomPosn[2]*k = 0)
        firstHydrogenRotationAxis = np.array([1, 1, -(axisToBond[0] + axisToBond[1])/axisToBond[2]])
        firstHydrogenRotationAxis /= np.linalg.norm(firstHydrogenRotationAxis)

        [a, b, c] = currentAtomPosn
        [u, v, w] = firstHydrogenRotationAxis
        # First hydrogen
        theta = 70.5 * np.pi/180.0
        newPosition = np.array([(a * (v**2 + w**2) - u * ((b * v) + (c * w) - (u * x) - (v * y) - (w * z))) * (1 - np.cos(theta)) + (x * np.cos(theta)) + ((-(c * v) + (b * w) - (w * y) + (v * z)) * np.sin(theta)),
                           (b * (u**2 + w**2) - v * ((a * u) + (c * w) - (u * x) - (v * y) - (w * z))) * (1 - np.cos(theta)) + (y * np.cos(theta)) + (((c * u)  - (a * w) + (w * x) - (u * z)) * np.sin(theta)),
                           (c * (u**2 + v**2) - w * ((a * u) + (b * v) - (u * x) - (v * y) - (w * z))) * (1 - np.cos(theta)) + (z * np.cos(theta)) + ((-(b * u) + (a * v) - (v * x) + (u * y)) * np.sin(theta))])
        hydrogenPositions.append(newPosition)
        # Second and third hydrogens
        # Rotate these from the newPosition +/-120 degrees around the vector axisToBond from the position currentAtomPosn - axisToBond
        [x, y, z] = newPosition
        [a, b, c] = currentAtomPosn + (np.cos(theta) * axisToBond)
        [u, v, w] = ((np.cos(theta) * axisToBond) / np.linalg.norm(np.cos(theta) * axisToBond))
        for theta in [120 * (np.pi / 180.0), -120 * (np.pi / 180.0)]:
            newHydrogenPosition = np.array([(a * (v**2 + w**2) - u * ((b * v) + (c * w) - (u * x) - (v * y) - (w * z))) * (1 - np.cos(theta)) + (x * np.cos(theta)) + ((-(c * v) + (b * w) - (w * y) + (v * z)) * np.sin(theta)),
                               (b * (u**2 + w**2) - v * ((a * u) + (c * w) - (u * x) - (v * y) - (w * z))) * (1 - np.cos(theta)) + (y * np.cos(theta)) + (((c * u)  - (a * w) + (w * x) - (u * z)) * np.sin(theta)),
                               (c * (u**2 + v**2) - w * ((a * u) + (b * v) - (u * x) - (v * y) - (w * z))) * (1 - np.cos(theta)) + (z * np.cos(theta)) + ((-(b * u) + (a * v) - (v * x) + (u * y)) * np.sin(theta))])
            hydrogenPositions.append(newHydrogenPosition)
    return hydrogenPositions


def loadMorphologyXMLETree(xmlPath, sigma=1.0):
    atomProps3DFloat = ['position']
    atomProps3DInt = ['image']
    atomPropsInt = ['body']
    atomPropsFloat = ['mass', 'diameter', 'charge']
    atomPropsStr = ['type']
    constraintProps = ['bond', 'angle', 'dihedral', 'improper']
    atomDictionary = {}
    with open(xmlPath, 'r') as xmlFileName:
        xmlFile = ET.parse(xmlFileName)
    morphologyConfig = xmlFile.getroot()[-1]
    for axis, systemDim in morphologyConfig.find('box').attrib.items():
        atomDictionary[axis] = float(systemDim)
    for key in atomPropsInt:
        atomDictionary[key] = list(map(int, morphologyConfig.find(key).text.split('\n')[1:-1]))
    for key in atomPropsFloat:
        atomDictionary[key] = list(map(float, morphologyConfig.find(key).text.split('\n')[1:-1]))
    for key in atomPropsStr:
        atomDictionary[key] = morphologyConfig.find(key).text.split('\n')[1:-1]
    for key in atomProps3DInt:
        atomDictionary[key] = [list(map(int, x.split(' '))) for x in morphologyConfig.find(key).text.split('\n')[1:-1]]
    for key in atomProps3DFloat:
        atomDictionary[key] = [list(np.array(list(map(float, x.split(' ')))) * sigma) for x in morphologyConfig.find(key).text.split('\n')[1:-1]]
    for key in constraintProps:
        atomDictionary[key] = [[x.split(' ')[0]] + list(map(int, x.split(' ')[1:])) for x in morphologyConfig.find(key).text.split('\n')[1:-1]]
    return atomDictionary


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
    AtomDictionary = {'position': [], 'image': [], 'mass': [], 'diameter': [], 'type': [], 'body': [], 'bond': [], 'angle': [], 'dihedral': [], 'improper': [], 'charge': []}
    record = False
    with open(xmlPath, 'r') as xmlFile:
        xmlData = xmlFile.readlines()
        for line in xmlData:
            if ('</' in line) or ('<!--' in line):
                record = False
            elif ('<configuration' in line) or ('<box' in line):
                # Get configuration data from this line (timestep, natoms etc)
                splitLine = line.split(' ')
                for i in range(1, len(splitLine)):
                    equalsLoc = findIndex(splitLine[i], '=')
                    if (equalsLoc is None) or ('units' in splitLine[i]):
                        # Skip any elements without equals or to do with the distance units
                        continue
                    quotationLoc = findIndex(splitLine[i], '"')
                    if ('.' in splitLine[i][quotationLoc[0] + 1:quotationLoc[1]]):
                        # Catch float in the value (excludes the = and quotation marks)
                        AtomDictionary[splitLine[i][:equalsLoc[0]].lower()] = float(splitLine[i][quotationLoc[0] + 1:quotationLoc[1]])
                    else:
                        AtomDictionary[splitLine[i][:equalsLoc[0]].lower()] = int(splitLine[i][quotationLoc[0] + 1:quotationLoc[1]])
            elif ('<position' in line):
                record = True
                recordType = 'position'
                continue
            elif ('<image' in line):
                record = True
                recordType = 'image'
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
            elif ('<bond' in line) and ('_coeff' not in line):
                record = True
                recordType = 'bond'
                continue
            elif ('<angle' in line) and ('_coeff' not in line):
                record = True
                recordType = 'angle'
                continue
            elif ('<dihedral' in line) and ('_coeff' not in line):
                record = True
                recordType = 'dihedral'
                continue
            elif ('<improper' in line) and ('_coeff' not in line):
                record = True
                recordType = 'improper'
                continue
            elif ('<charge' in line):
                record = True
                recordType = 'charge'
                continue
            # Now we know what the variable is, append it to the dictionary data
            if (record is True):
                # Mbuild outputs properties that are split by \t, so do a bit of jiggery pokery to allow us to interpret both
                splitLine = ' '.join(line.split('\t')).split(' ')
                # Remove the "\n"
                splitLine[-1] = splitLine[-1][:-1]
                if (recordType == 'position'):
                    # NOTE: VELOCITIES ARE NOT NORMALISED IN THE MORPHOLOGY FILE...DO THEY NEED TO BE SCALED BY SIGMA OR NOT? CURRENTLY THEY ARE.
                    # Write to dictionary as floats scaled by sigma
                    if (len(splitLine) == 1):
                        AtomDictionary[recordType].append(float(splitLine[0]))
                        continue
                    for i in range(len(splitLine)):
                        splitLine[i] = float(splitLine[i])
                    AtomDictionary[recordType].append(splitLine)
                elif (recordType == 'mass') or (recordType == 'diameter') or (recordType == 'charge'):
                    # Write to dictionary as floats
                    if (len(splitLine) == 1):
                        AtomDictionary[recordType].append(float(splitLine[0]))
                        continue
                    for i in range(len(splitLine)):
                        splitLine[i] = float(splitLine[i])
                    AtomDictionary[recordType].append(splitLine)
                elif (recordType == 'image') or (recordType == 'body'):
                    # Write to dictionary as int
                    if (len(splitLine) == 1):
                        AtomDictionary[recordType].append(int(splitLine[0]))
                        continue
                    for i in range(len(splitLine)):
                        splitLine[i] = int(splitLine[i])
                    AtomDictionary[recordType].append(splitLine)
                elif (recordType == 'type'):
                    # Write to dictionary as str
                    AtomDictionary[recordType].append(str(splitLine[0]))
                else:
                    #  (recordType == 'bond') or (recordType == 'angle') or (recordType == 'dihedral') or (recordType == 'improper')
                    # Write to dictionary as combination
                    splitLine[0] = str(splitLine[0])
                    for i in range(1, len(splitLine)):
                        splitLine[i] = int(splitLine[i])
                    AtomDictionary[recordType].append(splitLine)
    if sigma != 1.0:
        AtomDictionary = scale(AtomDictionary, sigma)
    return AtomDictionary


def loadFFXML(xmlPath, mapping = False):
    FFDict = {'lj':[], 'dpd':[], 'bond':[], 'angle':[], 'dihedral':[], 'improper':[]}
    with open(xmlPath, 'r') as xmlFile:
        xmlData = xmlFile.readlines()
        for line in xmlData:
            if ('</' in line):
                record = False
            elif ('<lj' in line):
                record = True
                recordType = 'lj'
                continue
            elif ('<dpd' in line):
                record = True
                recordType = 'dpd'
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
            # Now we know what the variable is, append it to the dictionary data
            if (record is True):
                # Write to dictionary as combination
                splitLine = line.split(' ')
                # Remove the "\n"
                splitLine[-1] = splitLine[-1][:-1]
                splitLine[0] = str(splitLine[0])
                for i in range(1, len(splitLine)):
                    splitLine[i] = float(splitLine[i])
                FFDict[recordType].append(splitLine)
    # Now remap the names of the constraints if any mappings have been specified
    if mapping is not False:
        for constraintType in list(FFDict.keys()):
            for index, constraint in enumerate(FFDict[constraintType]):
                # Split the constraint name up based on each atom type
                constraintName = copy.deepcopy(constraint[0].split('-'))
                # Remap each atom
                for atomLoc, atomType in enumerate(constraintName):
                    constraintName[atomLoc] = mapping[atomType]
                # Apply the mapping to the FFDict
                FFDict[constraintType][index][0] = '-'.join(constraintName)
    return FFDict


def checkConstraintNames(AAMorphologyDict):
    # A function that renames the constraints based on the atom types given in the dictionary
    constraintTypes = ['bond', 'angle', 'dihedral', 'improper']
    for constraintType in constraintTypes:
        for constraintID, constraint in enumerate(AAMorphologyDict[constraintType]):
            newConstraintName = ""
            # Iterate through the atomIDs and update the constraint name based on the types
            for atomID in constraint[1:]:
                newConstraintName += AAMorphologyDict['type'][atomID]
                newConstraintName += "-"
            # Update the dict if the name has changed
            if (constraint[0] != newConstraintName[:-1]):
                AAMorphologyDict[constraintType][constraintID][0] = newConstraintName[:-1]
    return AAMorphologyDict


def writeMorphologyXMLETree(inputDictionary, outputFile):
    print("\n \n THIS DOES NOT SUPPORT TILT FACTORS AT ALL!!!!!!!!!!! \n \n")
    print("Checking wrapped positions before writing XML...")
    inputDictionary = checkWrappedPositions(inputDictionary)
    systemProps = ['box']
    atomProps3D = ['position', 'image']
    atomProps = ['mass', 'diameter', 'type', 'body', 'charge']
    constraintProps = ['bond', 'angle', 'dihedral', 'improper']
    root = ET.Element('hoomd_xml', version="1.5")
    root.text = '\n'
    config = ET.Element('configuration', time_step=str(inputDictionary['time_step']), dimensions="3", natoms=str(inputDictionary['natoms']))
    config.text = '\n'
    config.tail = '\n'
    for element in systemProps + atomProps3D + atomProps + constraintProps:
        ET.SubElement(config, element)
        config[-1].text = '\n'
        config[-1].tail = '\n'
    for axis in ['lx', 'ly', 'lz']:
        config.find('box').attrib[axis] = str(inputDictionary[axis])
    for axis in ['xy', 'xz', 'yz']:
        config.find('box').attrib[axis] = str(0)
    config.find('box').text = ""
    config.attrib['natoms'] = str(inputDictionary['natoms'])
    for atomID, atomType in enumerate(inputDictionary['type']):
        for atomProp3D in atomProps3D:
            config.find(atomProp3D).text += ' '.join([str(x) for x in inputDictionary[atomProp3D][atomID]]) + '\n'
            config.find(atomProp3D).attrib['num'] = str(len(inputDictionary[atomProp3D]))
        for atomProp in atomProps:
            config.find(atomProp).text += str(inputDictionary[atomProp][atomID]) + '\n'
            config.find(atomProp).attrib['num'] = str(len(inputDictionary[atomProp]))
    for constraintType in constraintProps:
        for constraintID, constraint in enumerate(inputDictionary[constraintType]):
            config.find(constraintType).text += ' '.join([str(x) for x in inputDictionary[constraintType][constraintID]]) + '\n'
        config.find(constraintType).attrib['num'] = str(len(inputDictionary[constraintType]))
    root.insert(0, config)
    tree = ET.ElementTree(root)
    tree.write(outputFile, xml_declaration=True, encoding='UTF-8')
    print("XML file written to", str(outputFile) + "!")


def writeMorphologyXML(inputDictionary, outputFile, sigma = 1.0, checkWrappedPosns = True):
    tilt_factors = ["xy", "xz", "yz"]
    # Firstly, scale everything by the inverse of the provided sigma value
    if sigma != 1.0:
        inputDictionary = scale(inputDictionary, 1.0 / sigma)
    # Now need to check the positions of the atoms to ensure that everything is correctly contained inside the box
    if checkWrappedPosns is True:
        if any([tilt_factor in inputDictionary.keys() for tilt_factor in tilt_factors]) and\
                any([inputDictionary[tilt_factor] != 0 for tilt_factor in tilt_factors]):
            print("Can't check atom wrapping for cells with a non-zero tilt factor")
        else:
            print("Checking wrapped positions before writing XML...")
            inputDictionary = checkWrappedPositions(inputDictionary)
    # inputDictionary['position'], inputDictionary['image'] = pbc.shift_pbc(inputDictionary['position'], [inputDictionary['lx'], inputDictionary['ly'], inputDictionary['lz']])
    # print inputDictionary['image'][:20]
    # raw_input('HALT')
    # Add Boiler Plate first
    linesToWrite = ['<?xml version="1.0" encoding="UTF-8"?>\n', '<hoomd_xml version="1.4">\n', '<configuration time_step="0" dimensions="3" natoms="' + str(inputDictionary['natoms']) + '" >\n', '<box lx="' + str(inputDictionary['lx']) + '" ly="' + str(inputDictionary['ly']) + '" lz="' + str(inputDictionary['lz'])]
    if all([tilt_factor in inputDictionary.keys() for tilt_factor in tilt_factors]):
        linesToWrite[-1] += '" xy="' + str(inputDictionary['xy']) + '" xz="' + str(inputDictionary['xz']) + '" yz="' + str(inputDictionary['yz']) + '" />\n'
    else:
        linesToWrite[-1] += '" />\n'
    # Position
    linesToWrite.append('<position num="' + str(inputDictionary['natoms']) + '">\n')
    for positionData in inputDictionary['position']:
        linesToWrite.append(" ".join(str(coord) for coord in positionData) + '\n')
    linesToWrite.append('</position>\n')
    # Image
    linesToWrite.append('<image num="' + str(inputDictionary['natoms']) + '">\n')
    for imageData in inputDictionary['image']:
        linesToWrite.append(" ".join(str(coord) for coord in imageData) + '\n')
    linesToWrite.append('</image>\n')
    # Mass
    linesToWrite.append('<mass num="' + str(inputDictionary['natoms']) + '">\n')
    for massData in inputDictionary['mass']:
        linesToWrite.append(str(massData) + '\n')
    linesToWrite.append('</mass>\n')
    # Diameter
    linesToWrite.append('<diameter num="' + str(inputDictionary['natoms']) + '">\n')
    for diameterData in inputDictionary['diameter']:
        linesToWrite.append(str(diameterData) + '\n')
    linesToWrite.append('</diameter>\n')
    # Type
    linesToWrite.append('<type num="' + str(inputDictionary['natoms']) + '">\n')
    for typeData in inputDictionary['type']:
        linesToWrite.append(str(typeData) + '\n')
    linesToWrite.append('</type>\n')
    # Body
    linesToWrite.append('<body num="' + str(inputDictionary['natoms']) + '">\n')
    for bodyData in inputDictionary['body']:
        linesToWrite.append(str(bodyData) + '\n')
    linesToWrite.append('</body>\n')
    # Bond
    linesToWrite.append('<bond num="' + str(len(inputDictionary['bond'])) + '">\n')
    for bondData in inputDictionary['bond']:
        linesToWrite.append(" ".join(str(coord) for coord in bondData) + '\n')
    linesToWrite.append('</bond>\n')
    # Angle
    linesToWrite.append('<angle num="' + str(len(inputDictionary['angle'])) + '">\n')
    for angleData in inputDictionary['angle']:
        linesToWrite.append(" ".join(str(coord) for coord in angleData) + '\n')
    linesToWrite.append('</angle>\n')
    # Dihedral
    linesToWrite.append('<dihedral num="' + str(len(inputDictionary['dihedral'])) + '">\n')
    for dihedralData in inputDictionary['dihedral']:
        linesToWrite.append(" ".join(str(coord) for coord in dihedralData) + '\n')
    linesToWrite.append('</dihedral>\n')
    # Improper
    linesToWrite.append('<improper num="' + str(len(inputDictionary['improper'])) + '">\n')
    for improperData in inputDictionary['improper']:
        linesToWrite.append(" ".join(str(coord) for coord in improperData) + '\n')
    linesToWrite.append('</improper>\n')
    # Charge
    linesToWrite.append('<charge num="' + str(inputDictionary['natoms']) + '">\n')
    for chargeData in inputDictionary['charge']:
        linesToWrite.append(str(chargeData) + '\n')
    linesToWrite.append('</charge>\n')
    linesToWrite.append('</configuration>\n')
    linesToWrite.append('</hoomd_xml>\n')
    with open(outputFile, 'w+') as xmlFile:
        xmlFile.writelines(linesToWrite)
    print("XML file written to", str(outputFile) + "!")


def writeXYZFile(inputDict, outputFile):
    '''This function takes an input dictionary and converts it to an XYZ for use in DFT calculations'''
    # First line is atom numbers, second line is boiler plate
    rowsToWrite = [str(inputDict['natoms']) + '\n', 'XYZ file generated from XML using helperFunctions.XMLToXYZ\n']
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
        lineToWrite = atomType + atomX + atomY + atomZ + '\n'
        rowsToWrite.append(lineToWrite)
    with open(outputFile, 'w+') as xyzFile:
        xyzFile.writelines(rowsToWrite)
    print("XYZ data written to", str(outputFile) + ".")


def incrementAtomIDs(originalInputDictionary, ghostDictionary, increment, modifyGhostDictionary=False):
    inputDictionary = copy.deepcopy(originalInputDictionary)
    constraintTypes = ['bond', 'angle', 'dihedral', 'improper']
    for constraintType in constraintTypes:
        for constraintNo, constraint in enumerate(inputDictionary[constraintType]):
            inputDictionary[constraintType][constraintNo][1:] = [x + increment for x in inputDictionary[constraintType][constraintNo][1:]]
    if modifyGhostDictionary is True:
        for bondNo, bond in enumerate(ghostDictionary['bond']):
            if str(bond[1])[0] == '_':
                ghostDictionary['bond'][bondNo][1] = int(bond[1][1:]) + increment
            if str(bond[2])[0] == '_':
                ghostDictionary['bond'][bondNo][2] = int(bond[2][1:]) + increment
    return inputDictionary, ghostDictionary


def scale(inputDictionary, scaleFactor):
    for ID, position in enumerate(inputDictionary['position']):
        # if ID == 24104:
        #     print "Initial Position =", inputDictionary['position'][ID], inputDictionary['image'][ID]
        inputDictionary['position'][ID] = list(scaleFactor * np.array(position))
        # if ID == 24104:
        #     print "Scaled Position =", inputDictionary['position'][ID], inputDictionary['image'][ID]
    for element in ['lx', 'ly', 'lz']:
        if element in inputDictionary:
            inputDictionary[element] *= scaleFactor
    return inputDictionary


def rotate(inputDictionary, theta, rotateAroundPoint = [0, 0, 0], rotateAroundAxis = [0, 0, 1]):
    inputDictionary = addUnwrappedPositions(inputDictionary)
    rotateAroundAxis = list(np.array(rotateAroundAxis) / np.linalg.norm(rotateAroundAxis))
    # Rotation matrix calculations from: http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
    # The array that describes the 3D rotation of (x, y, z) around the point (a, b, c) through
    # the unit axis <u, v, w> by the angle theta is given by:
    # [ (a(v^2 + w^2) - u(bv + cw - ux - vy - wz))(1 - cos(theta)) + x*cos(theta) + (-cv + bw - wy + vz)sin(theta),
    #   (b(u^2 + w^2) - v(au + cw - ux - vy - wz))(1 - cos(theta)) + y*cos(theta) + (cu - aw + wx - uz)sin(theta),
    #   (c(u^2 + v^2) - w(au + bv - ux - vy - wz))(1 - cos(theta)) + z*cos(theta) + (-bu + av - vx + uy)sin(theta) ]
    # DEFAULT BEHAVIOUR: Rotate the entire dictionary by theta around the z-axis centred at the origin
    for AAID, [x, y, z] in enumerate(inputDictionary['unwrapped_position']):
        [a, b, c] = rotateAroundPoint
        [u, v, w] = rotateAroundAxis
        newPosition = np.array([(a * (v**2 + w**2) - u * ((b * v) + (c * w) - (u * x) - (v * y) - (w * z))) * (1 - np.cos(theta)) + (x * np.cos(theta)) + ((-(c * v) + (b * w) - (w * y) + (v * z)) * np.sin(theta)),
                       (b * (u**2 + w**2) - v * ((a * u) + (c * w) - (u * x) - (v * y) - (w * z))) * (1 - np.cos(theta)) + (y * np.cos(theta)) + (((c * u)  - (a * w) + (w * x) - (u * z)) * np.sin(theta)),
                       (c * (u**2 + v**2) - w * ((a * u) + (b * v) - (u * x) - (v * y) - (w * z))) * (1 - np.cos(theta)) + (z * np.cos(theta)) + ((-(b * u) + (a * v) - (v * x) + (u * y)) * np.sin(theta))])
        inputDictionary['unwrapped_position'][AAID] = list(newPosition)
    # All the images are probably messed up now, so fix that
    inputDictionary = replaceWrappedPositions(inputDictionary)
    return inputDictionary


def centre(inputDictionary, centreOfMass):
    COM = np.array(centreOfMass)
    for index, position in enumerate(inputDictionary['position']):
        inputDictionary['position'][index] = list(position - COM)
    return inputDictionary


def checkWrappedPositions(inputDictionary):
    atomPositions = np.array(inputDictionary['position'])
    atomImages = np.array(inputDictionary['image'])
    xhi = inputDictionary['lx'] / 2.0
    xlo = -inputDictionary['lx'] / 2.0
    yhi = inputDictionary['ly'] / 2.0
    ylo = -inputDictionary['ly'] / 2.0
    zhi = inputDictionary['lz'] / 2.0
    zlo = -inputDictionary['lz'] / 2.0
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


def getCPUCores():
    # Determine the number of available processors, either by querying the SLURM_NPROCS environment variable, or by using multiprocessing to count the number of visible CPUs.
    try:
        procIDs = list(np.arange(int(os.environ.get('SLURM_NPROCS'))))
    except (AttributeError, TypeError):
        # Was not loaded using SLURM, so use all physical processors
        procIDs = list(np.arange(mp.cpu_count()))
    return procIDs


def writeToFile(logFile, stringList, mode='logFile'):
    if mode == 'outputFile':
        openAs = 'w+'
    else:
        openAs = 'a+'
    if logFile == 'stdout':
        for line in stringList:
            sys.stdout.writelines(line + '\n')
    else:
        with open(logFile, openAs) as logWrite:
            for line in stringList:
                logWrite.writelines(line + '\n')


def loadPickle(pickleLocation):
    print("Loading Pickle from", str(pickleLocation) + "...")
    try:
        with open(pickleLocation, 'rb') as pickleFile:
                objects = pickle.load(pickleFile)
    except UnicodeDecodeError:  # Python 2/3 fix
        print("Old pickle! Loading it using Python 2...")
        with open(pickleLocation, 'rb') as pickleFile:
            objects = pickle.load(pickleFile, encoding='latin1')
    print("Pickle loaded successfully!")
    return objects[0], objects[1], objects[2], objects[3], objects[4]


def writePickle(toPickle, pickleFileName):
    print("Writing pickle file...")
    with open(pickleFileName, 'wb+') as pickleFile:
        pickle.dump(toPickle, pickleFile)
    print("Pickle file written to", pickleFileName)


def obtainBondedList(bondList):
    # Create a lookup table `neighbour list' for all connected atoms called {bondedAtoms}
    bondedAtoms = {}
    for bond in bondList:
        if bond[1] not in bondedAtoms:
            bondedAtoms[bond[1]] = [bond[2]]
        else:
            bondedAtoms[bond[1]].append(bond[2])
        if bond[2] not in bondedAtoms:
            bondedAtoms[bond[2]] = [bond[1]]
        else:
            bondedAtoms[bond[2]].append(bond[1])
    return bondedAtoms


def convertStringToInt(x):
    for i in range(len(x)):
        try:
            return int(x[i:])
        except:
            continue
    return 99999


def fixImages(originalMorphology):
    def checkBonds(morphology, bondDict):
        periodicBonds = []
        for bond in morphology['bond']:
            posn1 = np.array(morphology['position'][bond[1]]) + (np.array(morphology['image'][bond[1]]) * np.array([morphology['lx'], morphology['ly'], morphology['lz']]))
            posn2 = np.array(morphology['position'][bond[2]]) + (np.array(morphology['image'][bond[2]]) * np.array([morphology['lx'], morphology['ly'], morphology['lz']]))
            separation = calculateSeparation(posn1, posn2)
            if separation >= morphology['lx'] / 2.0:
                print("Periodic bond found:", bond, "because separation =", separation, ">=", morphology['lx'] / 2.0)
                morphology = moveBondedAtoms(bond[1], morphology, bondDict)
        return morphology

    def zeroOutImages(morphology):
        for atomID, image in enumerate(morphology['image']):
            if image != [0, 0, 0]:
                morphology['image'][atomID] = [0, 0, 0]
        return morphology

    def getBondDict(morphology):
        bondDict = {atomID: [] for atomID, atomType in enumerate(morphology['type'])}
        for bond in morphology['bond']:
            #if bond[1] < bond[2]:
            bondDict[bond[1]].append(bond[2])
            #else:
            bondDict[bond[2]].append(bond[1])
        return bondDict

    def moveBondedAtoms(centralAtom, morphology, bondDict):
        for bondedAtom in bondDict[centralAtom]:
            atom1Posn = morphology['position'][centralAtom]
            atom2Posn = morphology['position'][bondedAtom]
            #print("atom1:", centralAtom, "posn =", atom1Posn, "; atom2:", bondedAtom, "posn =", atom2Posn)
            sepVec = np.array(atom1Posn) - np.array(atom2Posn)
            moved = False
            for axis, value in enumerate(sepVec):
                if value > morphology['lx'] / 2.0:
                    morphology['position'][bondedAtom][axis] += morphology['lx']
                    moved = True
                if value < -morphology['lx'] / 2.0:
                    morphology['position'][bondedAtom][axis] -= morphology['lx']
                    moved = True
            if moved:
                #print("Moved", bondedAtom, "to same box as", centralAtom)
                #print("New Positions: atom1 =", morphology['position'][centralAtom], "atom2 =", morphology['position'][bondedAtom])
                morphology = moveBondedAtoms(bondedAtom, morphology, bondDict)
            else:
                #print("Move was unnecessary")
                pass
        return morphology

    zeroedMorphology = zeroOutImages(originalMorphology)
    bondDict = getBondDict(zeroedMorphology)
    fixedMorphology = checkBonds(zeroedMorphology, bondDict)
    return fixedMorphology


# ---============================---
# ---=== KMC HELPER FUNCTIONS ===---
# ---============================---
def calculateCarrierHopRate(lambdaij, Tij, deltaEij, prefactor, temp, useVRH=False, rij=0.0, VRHPrefactor=1.0, boltzPen=False):
    # Based on the input parameters, can make this the semiclassical Marcus Hopping Rate Equation, or a more generic Miller Abrahams-based hop
    # Firstly, to prevent divide-by-zero errors:
    if (Tij == 0.0):
        return 0
    # Regardless of hopping type, sort out the prefactor first:
    kij = prefactor * ((2 * np.pi) / hbar) * (Tij ** 2) * np.sqrt(1.0 / (4 * lambdaij * np.pi * kB * temp))
    # VRH?
    if useVRH is True:
        kij *= np.exp(-(VRHPrefactor * rij))
    # Simple Boltzmann energy penalty?
    if boltzPen is True:
        # Only apply the penalty if deltaEij is positive
        if deltaEij > 0.0:
            kij *= np.exp(-(deltaEij / (kB * temp)))
        # Otherwise, kij *= 1
    else:
        kij *= np.exp(-((deltaEij + lambdaij)**2) / (4 * lambdaij * kB * temp))
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


def calculateMillerAbrahamsHopRate(prefactor, separation, radius, deltaEij, T):
    kij = prefactor * np.exp(-2 * separation/radius)
    if deltaEij > 0:
        kij *= np.exp(-deltaEij / (kB * T))
    return kij


def determineEventTau(rate, eventType='None', slowestEvent=None, fastestEvent=None, maximumAttempts=None):
    # Use the KMC algorithm to determine the wait time to this hop
    if rate != 0:
        counter = 0
        while True:
            if maximumAttempts is not None:
                # Write an error if we've hit the maximum number of attempts
                if counter == maximumAttempts:
                    if 'hop' in eventType:
                        return None
                    else:
                        if 'injection' not in eventType:
                            helperFunctions.writeToFile(logFile, ["Attempted " + str(maximumAttempts) + " times to obtain a '" +
                                                                  str(eventType) + "'-type event timescale within the tolerances: " +
                                                                  str(fastestEvent) + " <= tau < " +
                                                                  str(slowestEvent) + " with the given rate " +
                                                                  str(rate) + " all without success.",
                                                                  "Permitting the event anyway with the next random number."])

            x = np.random.random()
            # Ensure that we don't get exactly 0.0 or 1.0, which would break our logarithm
            if (x == 0.0) or (x == 1.0):
                continue
            tau = - np.log(x) / rate
            if (fastestEvent is not None) and (slowestEvent is not None) and (maximumAttempts is not None):
                if ((tau > fastestEvent) and (tau < slowestEvent)) or (counter == maximumAttempts):
                    break
                else:
                    counter += 1
                    continue
            break
    else:
        # If rate == 0, then make the hopping time extremely long
        tau = 1E99
    return tau
