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


def calculateHydrogenPositions(morphologyDict, atomsToAddTo):
    '''This function calculates the position of the hydrogen based
    on the positions of the other bonded species'''
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
    print("THIS FUNCTION IS NOT YET GENERAL PURPOSE! IT WORKS ONLY FOR AROMATIC SMALL MOLECULES WHERE THE CARBONS HAVE ONE SINGLE AND ONE DOUBLE BOND EXACTLY")
    for inputFile in os.listdir('./'):
        if ('_AA.xml' in inputFile) or ('.xml' not in inputFile):
            continue
        morphologyDict = helperFunctions.loadMorphologyXML(inputFile, sigma = 3.8)
        morphologyDict = helperFunctions.addUnwrappedPositions(morphologyDict)
        atomsToAddTo = hydrogenIndices(morphologyDict)
        hydrogenPositions = calculateHydrogenPositions(morphologyDict, atomsToAddTo)
        morphologyDict = addHydrogensToMorph(morphologyDict, hydrogenPositions)
        morphologyDict = helperFunctions.addWrappedPositions(morphologyDict)
        helperFunctions.writeMorphologyXML(morphologyDict, inputFile.replace('.xml','_AA.xml'))
