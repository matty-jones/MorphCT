import numpy as np
import os
import sys
sys.path.append('../../code')
import helperFunctions


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
    for atomID, atomType in enumerate(morphologyDict['type']):
        # Skip if we don't have to add hydrogens to the current atom's type
        if atomType not in hydrogensToAdd.keys():
            continue
        for bondDefinition in hydrogensToAdd[atomType]:
            # Skip if the current atom does not have the right number of bonds
            if numberOfBonds[atomID][0] != bondDefinition[0]:
                continue
            # Otherwise, we need to add hydrogensToAdd[atomType][1] hydrogens to this atom
            currentAtomPosn = morphologyDict['unwrapped_position'][atomID]
            # First get the vector to the average position of the bonded neighbours
            averagePositionOfBondedAtoms = np.array([0.0, 0.0, 0.0])
            for bondedAtom in numberOfBonds[atomID][1]:
                bondVector = np.array(morphologyDict['unwrapped_position'][bondedAtom]) - currentAtomPosn
                bondVector /= np.linalg.norm(bondVector)
                averagePositionOfBondedAtoms += bondVector
            [x, y, z] = currentAtomPosn + (-1.06 * (averagePositionOfBondedAtoms / np.linalg.norm(averagePositionOfBondedAtoms)))
            if bondDefinition[1] == 1:
                # Easy, this is the perylene code
                # Simply reverse the bonded vector and make it the hydrogen position at a distance of 1.06 angstroems
                hydrogenPositions.append([int(atomID), np.array([x, y, z])])
            # Initial position for all hydrogens
            elif bondDefinition[1] == 2:
                # As above (to get the right plane), but then rotated +(109.5/2) degrees and -(109.5/2) degrees around the bonding axis
                rotationAxis = np.array(morphologyDict['unwrapped_position'][numberOfBonds[atomID][1][0]]) - np.array(morphologyDict['unwrapped_position'][numberOfBonds[atomID][1][-1]])
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
                    hydrogenPositions.append([int(atomID), newPosition])
            elif bondDefinition[1] == 3:
                # As for one (to get the right side of the bonded atom), rotate the first one up by 70.5 (180 - 109.5) and then rotate around by 109.5 degrees for the other two
                # The first hydrogen can be rotated around any axis perpendicular to the only bond present
                axisToBond = currentAtomPosn - np.array(morphologyDict['unwrapped_position'][numberOfBonds[atomID][1][0]])
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
                hydrogenPositions.append([int(atomID), newPosition])
                # Second and third hydrogens
                # Rotate these from the newPosition +/-120 degrees around the vector axisToBond from the position currentAtomPosn - axisToBond
                [x, y, z] = newPosition
                [a, b, c] = currentAtomPosn + (np.cos(theta) * axisToBond)
                [u, v, w] = ((np.cos(theta) * axisToBond) / np.linalg.norm(np.cos(theta) * axisToBond))
                for theta in [120 * (np.pi / 180.0), -120 * (np.pi / 180.0)]:
                    newHydrogenPosition = np.array([(a * (v**2 + w**2) - u * ((b * v) + (c * w) - (u * x) - (v * y) - (w * z))) * (1 - np.cos(theta)) + (x * np.cos(theta)) + ((-(c * v) + (b * w) - (w * y) + (v * z)) * np.sin(theta)),
                                       (b * (u**2 + w**2) - v * ((a * u) + (c * w) - (u * x) - (v * y) - (w * z))) * (1 - np.cos(theta)) + (y * np.cos(theta)) + (((c * u)  - (a * w) + (w * x) - (u * z)) * np.sin(theta)),
                                       (c * (u**2 + v**2) - w * ((a * u) + (b * v) - (u * x) - (v * y) - (w * z))) * (1 - np.cos(theta)) + (z * np.cos(theta)) + ((-(b * u) + (a * v) - (v * x) + (u * y)) * np.sin(theta))])
                    hydrogenPositions.append([int(atomID), newHydrogenPosition])
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
        morphologyDict['mass'].append(1.00794)
        morphologyDict['charge'].append(0.0)
        morphologyDict['diameter'].append(0.53)
        morphologyDict['body'].append(morphologyDict['body'][hydrogenAtom[0]])
    return morphologyDict


if __name__ == "__main__":
    # This dictionary has keys of the atom type, and values where the first element is the number of bonds required for us to add a hydrogen to the atom and the second element of the value defines how many hydrogens to add to said atom.
    # ---==== PCBM ====---
    #print("Using data for PCBM")
    #hydrogensToAdd = {'CHA':[[2, 1]], # If the atom type is CHA and has only 2 bonds, add 1 hydrogen
    #                  'CH2':[[2, 2]], # If the atom type is CH2 and has only 2 bonds, add 2 hydrogens
    #                  'CE':[[1, 3]]}  # If the atom type is CE and has only one bond, add 3 hydrogens
    #sigmaVal = 1.0
    # ---==== P3HT ====---
    #print("Using data for P3HT")
    #hydrogensToAdd = {'CA':[[2, 1]],
    #                  'CT':[[2, 2],[1, 3]]}
    #sigmaVal = 3.905
    # ---==== PERYLENE/PERYLOTHIOPHENE ====---
    #print("Using data for Perylene/Perylothiophene")
    #hydrogensToAdd = {'C':[[2, 1]]}
    #sigmaVal = 3.905
    ## ---==== BDT-TPD ====---
    print("Using data for BDT-TPD")
    hydrogensToAdd = {'CS':[[2, 1]],
                      'C!':[[2, 1]],
                      'CT':[[2, 2],[1, 3], [3, 1]],
                      'CP':[[2, 1]]}
    sigmaVal = 3.55

    print("THIS FUNCTION IS SET UP TO USE A DICTIONARY TO DEFINE HOW MANY HYDROGENS TO ADD TO BONDS OF A SPECIFIC TYPE WITH A CERTAIN NUMBER OF BONDS")
    print(hydrogensToAdd)
    print("IF THE ABOVE DICTIONARY DOESN'T LOOK RIGHT, PLEASE TERMINATE NOW AND IGNORE ANY OUTPUTS UNTIL THE DICTIONARY HAS BEEN RECTIFIED")
    print("Additionally, we're using a sigma value of", sigmaVal)
    for inputFile in os.listdir('./'):
        if ('_AA.xml' in inputFile) or ('.xml' not in inputFile):
            continue
        morphologyDict = helperFunctions.loadMorphologyXML(inputFile, sigma = sigmaVal)
        morphologyDict = helperFunctions.addUnwrappedPositions(morphologyDict)
        hydrogenPositions = calculateHydrogenPositions(morphologyDict, hydrogensToAdd)
        morphologyDict = addHydrogensToMorph(morphologyDict, hydrogenPositions)
        morphologyDict = helperFunctions.addWrappedPositions(morphologyDict)
        helperFunctions.writeMorphologyXML(morphologyDict, inputFile.replace('.xml','_AA.xml'), checkWrappedPosns=False)#, sigma = sigmaVal)
