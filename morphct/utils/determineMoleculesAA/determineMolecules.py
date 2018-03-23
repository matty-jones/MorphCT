import sys

sys.path.append('../../code')
import helperFunctions


def splitMolecules(inputDictionary):
    # Split the full morphology into individual molecules
    moleculeAAIDs = []
    moleculeLengths = []
    # Create a lookup table `neighbour list' for all connected atoms called {bondedAtoms}
    bondedAtoms = helperFunctions.obtainBondedList(inputDictionary['bond'])
    moleculeList = [i for i in range(len(inputDictionary['type']))]
    # Recursively add all atoms in the neighbour list to this molecule
    for molID in range(len(moleculeList)):
        moleculeList = updateMolecule(molID, moleculeList, bondedAtoms)
    # Create a dictionary of the molecule data
    moleculeData = {}
    for atomID in range(len(inputDictionary['type'])):
        if moleculeList[atomID] not in moleculeData:
            moleculeData[moleculeList[atomID]] = [atomID]
        else:
            moleculeData[moleculeList[atomID]].append(atomID)
    # Return the list of AAIDs and the lengths of the molecules
    for moleculeID in list(moleculeData.keys()):
        moleculeAAIDs.append(sorted(moleculeData[moleculeID]))
        moleculeLengths.append(len(moleculeData[moleculeID]))
    return moleculeAAIDs, moleculeLengths


def updateMolecule(atomID, moleculeList, bondedAtoms):
    # Recursively add all neighbours of atom number atomID to this molecule
    try:
        for bondedAtom in bondedAtoms[atomID]:
            # If the moleculeID of the bonded atom is larger than that of the current one,
            # update the bonded atom's ID to the current one's to put it in this molecule,
            # then iterate through all of the bonded atom's neighbours
            if moleculeList[bondedAtom] > moleculeList[atomID]:
                moleculeList[bondedAtom] = moleculeList[atomID]
                moleculeList = updateMolecule(bondedAtom, moleculeList, bondedAtoms)
            # If the moleculeID of the current atom is larger than that of the bonded one,
            # update the current atom's ID to the bonded one's to put it in this molecule,
            # then iterate through all of the current atom's neighbours
            elif moleculeList[bondedAtom] < moleculeList[atomID]:
                moleculeList[atomID] = moleculeList[bondedAtom]
                moleculeList = updateMolecule(atomID, moleculeList, bondedAtoms)
            # Else: both the current and the bonded atom are already known to be in this
            # molecule, so we don't have to do anything else.
    except KeyError:
        # This means that there are no bonded CG sites (i.e. it's a single molecule)
        pass
    return moleculeList


if __name__ == "__main__":
    pickleLoc = sys.argv[1]
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = helperFunctions.loadPickle(pickleLoc)
    # Create the new CGToAAIDMaster
    if CGToAAIDMaster is None:
        moleculeAAIDs, moleculeLengths = splitMolecules(AAMorphologyDict)
        print(moleculeAAIDs)
