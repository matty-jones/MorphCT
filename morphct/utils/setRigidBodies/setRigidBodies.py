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
    fileName = sys.argv[1]
    morphology = helperFunctions.loadMorphologyXML(fileName)
    ## Got to fix the awful mess that is the final XML from the P3HT_C60 interface.
    ## Do the P3HT first
    #rollingBodyNumber = 0
    #bodyMappings = {}
    #P3HTAtoms = ['C10', 'C2', 'C9', 'S1', 'C1']
    #PCBMAtoms = ['C11']
    #for atomID, bodyID in enumerate(morphology['body']):
    #    if (bodyID != -1):
    #        if morphology['type'][atomID] in P3HTAtoms:
    #            if bodyID not in bodyMappings:
    #                bodyMappings[bodyID] = rollingBodyNumber
    #                morphology['body'][atomID] = rollingBodyNumber
    #                rollingBodyNumber += 1
    #            else:
    #                morphology['body'][atomID] = bodyMappings[bodyID]
    ## Then do the PCBM
    #bodyMappings = {}
    #for atomID, bodyID in enumerate(morphology['body']):
    #    if (bodyID != -1):
    #        if morphology['type'][atomID] in PCBMAtoms:
    #            if bodyID not in bodyMappings :
    #                bodyMappings[bodyID] = rollingBodyNumber
    #                morphology['body'][atomID] = rollingBodyNumber
    #                rollingBodyNumber += 1
    #            else:
    #                morphology['body'][atomID] = bodyMappings[bodyID]

    ## Finally, write out the xml
    #helperFunctions.writeMorphologyXML(morphology, "bodyFix_" + fileName)

    ## ---=== Stuff for PCBM ===---
    #PCBMAtomIDs = [AAID for AAID, atomType in enumerate(morphology['type']) if atomType in ['FCA', 'FCT', 'O']]
    #atomsPerMolecule = 74
    ## ---======================---
    #numberOfMols = len(PCBMAtomIDs)//atomsPerMolecule
    #startIndex = min(PCBMAtomIDs)
    #endIndex = max(PCBMAtomIDs)

    ## Set the first molNo to be == the AAID to ensure that there are no conflicts possible
    #for atomID, rigidBody in enumerate(morphology['body']):
    #    # Go up to where the molecules we care about stop
    #    if atomID > endIndex:
    #        break
    #    # Start where the molecules we care about start
    #    if atomID >= startIndex:
    #        # Set the body to be equal to the molecule id
    #        morphology['body'][atomID] = startIndex + ((atomID - startIndex) // atomsPerMolecule)

    ## Now make the rigid bodies all sequential (doesn't really matter but keeps it neat)
    #sortedBodyList = sorted(list(set(morphology['body'])))
    ## Take into account flexible bodies by subtracting all IDs by one
    #if sortedBodyList[0] == -1:
    #    modifier = -1
    #else:
    #    modifier = 0
    #for atomID, rigidBody in enumerate(morphology['body']):
    #    morphology['body'][atomID] = sortedBodyList.index(rigidBody) + modifier


    # Set rigid bodies to be equal to the molID
    moleculeAAIDs, moleculeLengths = splitMolecules(morphology)
    AAIDToMolID = {}
    for index, moleculeAAIDList in enumerate(moleculeAAIDs):
        for AAID in moleculeAAIDList:
            AAIDToMolID[AAID] = index

    molIDs = []
    for AAID, atomType in enumerate(morphology['type']):
        # Get only the PCBM molecules:
        #if atomType == 'FCA':
        molIDs.append(AAIDToMolID[AAID])
    molIDs = list(set(molIDs))


    for AAID, molID in AAIDToMolID.items():
        if molID in molIDs:
            morphology['body'][AAID] = molID + len(morphology['body'])  # To prevent rigid body conflicts

    rigidBodies = sorted(list(set(morphology['body'])))
    try:
        rigidBodies.remove(-1)
    except ValueError:
        pass
    rigidBodyLookup = {-1: -1}
    for index, bodyNo in enumerate(rigidBodies):
        rigidBodyLookup[bodyNo] = index
    for atomID, rigidBodyID in enumerate(morphology['body']):
        morphology['body'][atomID] = rigidBodyLookup[rigidBodyID]


    # Finally, write out the xml
    helperFunctions.writeMorphologyXML(morphology, "bodyFix_" + fileName)
