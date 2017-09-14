import sys
sys.path.append('../../code')
import helperFunctions


if __name__ == "__main__":
    fileName = sys.argv[1]
    morphology = helperFunctions.loadMorphologyXML(fileName)
    # ---=== Stuff for PCBM ===---
    PCBMAtomIDs = [AAID for AAID, atomType in enumerate(morphology['type']) if atomType in ['FCA', 'FCT', 'O']]
    atomsPerMolecule = 74
    # ---======================---
    numberOfMols = len(PCBMAtomIDs)//atomsPerMolecule
    startIndex = min(PCBMAtomIDs)
    endIndex = max(PCBMAtomIDs)

    # Set the first molNo to be == the AAID to ensure that there are no conflicts possible
    for atomID, rigidBody in enumerate(morphology['body']):
        # Go up to where the molecules we care about stop
        if atomID > endIndex:
            break
        # Start where the molecules we care about start
        if atomID >= startIndex:
            # Set the body to be equal to the molecule id
            morphology['body'][atomID] = startIndex + ((atomID - startIndex) // atomsPerMolecule)

    # Now make the rigid bodies all sequential (doesn't really matter but keeps it neat)
    sortedBodyList = sorted(list(set(morphology['body'])))
    # Take into account flexible bodies by subtracting all IDs by one
    if sortedBodyList[0] == -1:
        modifier = -1
    else:
        modifier = 0
    for atomID, rigidBody in enumerate(morphology['body']):
        morphology['body'][atomID] = sortedBodyList.index(rigidBody) + modifier

    # Finally, write out the xml
    helperFunctions.writeMorphologyXML(morphology, "bodyFix_" + fileName)
