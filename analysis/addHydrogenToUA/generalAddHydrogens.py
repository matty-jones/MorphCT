import numpy as np
import os
import sys
sys.path.append('../../code')
import helperFunctions


def hydrogenIndices(morphologyDict):
    '''A function to determine the atomic indices to add hydrogens to.
    If a carbon atom has 2 bonded neighbours then add a hydrogen in a
    sensible position.'''
    bondList = morphologyDict['bond']
    atomsToAddTo = {}
    numberOfBonds = {}
    for bond in bondList:
        type1 = morphologyDict['type'][bond[1]]
        type2 = morphologyDict['type'][bond[2]]
        if (type1[0] == 'C') and (type2[0] == 'C'):
            if str(bond[1]) not in numberOfBonds:
                numberOfBonds[str(bond[1])] = []
            if str(bond[2]) not in numberOfBonds:
                numberOfBonds[str(bond[2])] = []
            numberOfBonds[str(bond[1])].append(morphologyDict['unwrapped_position'][bond[2]])
            numberOfBonds[str(bond[2])].append(morphologyDict['unwrapped_position'][bond[1]])
        elif (type1[0] == 'S'):
            if str(bond[2] not in numberOfBonds):
                numberOfBonds[str(bond[2])] = []
            numberOfBonds[str(bond[2])].append(morphologyDict['unwrapped_position'][bond[1]])
        elif (type2[0] == 'S'):
            if str(bond[1] not in numberOfBonds):
                numberOfBonds[str(bond[1])] = []
            numberOfBonds[str(bond[1])].append(morphologyDict['unwrapped_position'][bond[2]])
    for atomIndex, bondedPosns in numberOfBonds.items():
        if len(bondedPosns) == 2:
            atomsToAddTo[int(atomIndex)] = bondedPosns
    return atomsToAddTo


def calculateHydrogenPositions(morphologyDict, hydrogensToAdd):
    '''This function calculates the position of the hydrogen based
    on the number and positions of the other bonded species, and
    the number of hydrogens required to be added to the atom'''
    hydrogenPositions = []
    # First create a lookup table that describes exactly how many bonds
    # each atom has to neighbours
    numberOfBonds = {x:[0, []] for x in range(len(morphologyDict['type']))}  # Keys == atomID, values == [numberOfBondedSpecies, [bondedAtomIDs]]
    for bond in morphologyDict['bond']:
        atomID1 = bond[1]
        atomID2 = bond[2]
        if atomID2 not in numberOfBonds[atomID1][1]:
            numberOfBonds[atomID1][0] += 1
            numberOfBonds[atomID1][1].append(atomID2)
        if atomID1 not in numberOfBonds[atomID2][1]:
            numberOfBonds[atomID2][0] += 1
            numberOfBonds[atomID2][1].append(atomID1)
    for atomID, atomType in morphologyDict['type']:
        # Skip if we don't have to add hydrogens to the current atom's type
        if atomType not in hydrogensToAdd.keys():
            continue
        # Skip if the current atom does not have the right number of bonds
        if numberOfBonds[atomID][0] != hydrogensToAdd[atomType][0]:
            continue
        # Otherwise, we need to add hydrogensToAdd[atomType][1] hydrogens to this atom
        currentAtomPosn = morphologyDict['unwrapped_position'][atomID]
        # First get the vector to the average position of the bonded neighbours
        averagePositionOfBondedAtoms = np.array([0.0, 0.0, 0.0])
        for bondedAtom in numberOfBonds[atomID][1]:
            averagePositionOfBondedAtoms += np.array(morphologyDict['unwrapped_position'][bondedAtom]) - currentAtomPosn
        if hydrogensToAdd[atomType][1] == 1:
            # Easy, this is the perylene code
            # Simply reverse the bonded vector and make it the hydrogen position at a distance of 1.06 angstroems
            hydrogenPositions.append([int(atomID), currentAtomPosn + (-1.06 * (averagePositionOfBondedAtoms / np.linalg.norm(averagePositionOfBondedAtoms)))])
        elif hydrogensToAdd[atomType][1] == 2:
            # As above (to get the right plane), but then rotated +60 degrees and -60 degrees around the bonding axis
            rotationAxis = np.array(morphologyDict['unwrapped_position'][numberOfBonds[atomID][1][0]]) - np.array(morphologyDict['unwrapped_position'][numberOfBonds[atomID][1][-1]])
            rotationAxis /= np.linalg.norm(rotationAxis)
            # Rotation matrix calculations from: http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
            # The rotation matrix that describes the 3D rotation of (x, y, z) around the point (a, b, c) through
            # the unit axis <u, v, w> by the angle theta is given by:
            # [ [ (a(v^2 + w^2) - u(bv + cw - ux - vy - wz))(1 - cos(theta)) + x*cos(theta) + (-cv + bw - wy + vz)sin(theta) ],
            #   [ (b(u^2 + w^2) - v(au + cw - ux - vy - wz))(1 - cos(theta)) + y*cos(theta) + (cu - aw + wx - uz)sin(theta) ],
            #   [ (c(u^2 + v^2) - w(au + bv - ux - vy - wz))(1 - cos(theta)) + z*cos(theta) + (-bu + av - vx + uy)sin(theta) ] ]


            pass
        elif hydrogensToAdd[atomType][1] == 3:
            # As for one (to get the right side of the bonded atom) and then separated 109.5 degrees from each other
            pass


    exit()
    return hydrogenPositions





    hydrogenPositions = []
    for carbonAtom, bondedPositions in atomsToAddTo.items():
        carbonPosn = np.array(morphologyDict['unwrapped_position'][int(carbonAtom)])
        # Now find the position of the hydrogen.
        # First get the vector to the average position of the two bonded
        # neighbours...:
        averagePositionOfBondedAtoms = np.array([0.0, 0.0, 0.0])
        for position in bondedPositions:
            averagePositionOfBondedAtoms += np.array(position) - carbonPosn
        # Then reverse that vector and make it the hydrogen position at a
        # distance of 1.06 angstroems
        hydrogenPositions.append([int(carbonAtom), carbonPosn + (-1.06 * (averagePositionOfBondedAtoms / np.linalg.norm(averagePositionOfBondedAtoms)))])
    return hydrogenPositions


def addHydrogensToMorph(morphologyDict, hydrogenPositions):
    '''This function adds the hydrogen atoms into the morphology
    to be exported as the AA xml'''
    # HydrogenPositions are of the following format:
    # [CarbonIndex, PositionOfHydrogen]
    for hydrogenAtom in hydrogenPositions:
        morphologyDict['type'].append('H1')
        morphologyDict['unwrapped_position'].append(hydrogenAtom[1])
        morphologyDict['bond'].append(['C-H', hydrogenAtom[0], morphologyDict['natoms']])
        morphologyDict['natoms'] += 1
        otherProperties = ['mass', 'charge', 'body', 'diameter']
        for propName in otherProperties:
            morphologyDict[propName].append(morphologyDict[propName][hydrogenAtom[0]])
    return morphologyDict





if __name__ == "__main__":
    # This dictionary has keys of the atom type, and values where the first element is the number of bonds required for us to add a hydrogen to the atom and the second element of the value defines how many hydrogens to add to said atom.
    hydrogensToAdd = {'CHA':[2, 1], # If the atom type is CHA and has only 2 bonds, add 1 hydrogen
                      'CH2':[2, 2], # If the atom type is CH2 and has only 2 bonds, add 2 hydrogens
                      'CE':[1, 3]}  # If the atom type is CE and has only one bond, add 3 hydrogens

    print(hydrogensToAdd)
    print("THIS FUNCTION IS SET UP TO USE THE PREVIOUS DICTIONARY TO DEFINE HOW MANY HYDROGENS TO ADD TO BONDS OF A SPECIFIC TYPE WITH A CERTAIN NUMBER OF BONDS")
    print("IF THE ABOVE DICTIONARY DOESN'T LOOK RIGHT, PLEASE TERMINATE NOW AND IGNORE ANY OUTPUTS UNTIL THE DICTIONARY HAS BEEN RECTIFIED")
    for inputFile in os.listdir('./'):
        if ('_AA.xml' in inputFile) or ('.xml' not in inputFile):
            continue
        morphologyDict = helperFunctions.loadMorphologyXML(inputFile, sigma = 1.0)
        morphologyDict = helperFunctions.addUnwrappedPositions(morphologyDict)
        hydrogenPositions = calculateHydrogenPositions(morphologyDict, hydrogensToAdd)
        morphologyDict = addHydrogensToMorph(morphologyDict, hydrogenPositions)
        morphologyDict = helperFunctions.addWrappedPositions(morphologyDict)
        helperFunctions.writeMorphologyXML(morphologyDict, inputFile.replace('.xml','_AA.xml'))
