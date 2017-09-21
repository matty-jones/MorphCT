import sys
import numpy as np
sys.path.append('../../code')
import helperFunctions


def checkBonds(morphology, bondDict):
    periodicBonds = []
    for bond in morphology['bond']:
        posn1 = np.array(morphology['position'][bond[1]]) + (np.array(morphology['image'][bond[1]]) * np.array([morphology['lx'], morphology['ly'], morphology['lz']]))
        posn2 = np.array(morphology['position'][bond[2]]) + (np.array(morphology['image'][bond[2]]) * np.array([morphology['lx'], morphology['ly'], morphology['lz']]))
        separation = helperFunctions.calculateSeparation(posn1, posn2)
        if separation >= morphology['lx'] / 2.0:
            #print("Periodic bond found:", bond, "because separation =", separation, ">=", morphology['lx'] / 2.0)
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
    return morphology


if __name__ == "__main__":
    fileName = sys.argv[1]
    morphology = helperFunctions.loadMorphologyXML(fileName)
    morphology = zeroOutImages(morphology)
    bondDict = getBondDict(morphology)
    morphology = checkBonds(morphology, bondDict)
    helperFunctions.writeMorphologyXML(morphology, "fixed_" + fileName)
