from hoomd_script import *
import numpy as np
import copy
import sys
from morphct.code import helper_functions as hf


def obtainMoleculeDict(AAMorphologyDict, moleculeAAIDs):
    atomProps = ['position', 'image', 'mass', 'diameter', 'type', 'body', 'charge']
    constraintProps = ['bond', 'angle', 'dihedral', 'improper']
    systemProps = ['lx', 'ly', 'lz']
    # Synthesise an empty dictionary
    thisMoleculeDict = {}
    for key in atomProps + constraintProps:
        thisMoleculeDict[key] = []
    # Include the current system properties (needed to calculate unwrapped coords later)
    for key in systemProps:
        thisMoleculeDict[key] = AAMorphologyDict[key]
    # Create a dictionary to map the morphology AAIDs to the AAIDs within this molecule
    AAIDConversion = {}
    for newIndex, AAID in enumerate(moleculeAAIDs):
        AAIDConversion[str(AAID)] = newIndex
        for key in atomProps:
            thisMoleculeDict[key].append(AAMorphologyDict[key][AAID])
    # Iterate over all constraint types
    for constraintType in constraintProps:
        # For each constraint
        for constraint in AAMorphologyDict[constraintType]:
            # Check that all of the constraint IDs belong to this molecule
            # (Have to check all because ghost constraints are present in the AAMorphologyDict)
            if sum([x in constraint[1:] for x in moleculeAAIDs]) == len(constraint) - 1:
                newConstraint = copy.deepcopy(constraint)
                # Update all of the AAIDs to their new IDs within this molecule
                for location, atomID in enumerate(newConstraint):
                    if location == 0:
                        continue
                    newConstraint[location] = AAIDConversion[str(atomID)]
                # Now add them to the molecule dict
                thisMoleculeDict[constraintType].append(newConstraint)
    # Now need to make sure that everything is unwrapped
    thisMoleculeDict = hf.addUnwrappedPositions(thisMoleculeDict)
    # Find the (geometric, in case masses aren't specified) centre of the molecule
    positionArray = np.array(thisMoleculeDict['unwrapped_position'])
    centre = [np.average(positionArray[0, :]), np.average(positionArray[1, :]), np.average(positionArray[2, :])]
    # Also find the extent of the molecule so we can change lx, ly and lz
    boxDims = [np.max(positionArray[0, :]) - np.min(positionArray[0, :]), np.max(positionArray[1, :]) - np.min(positionArray[1, :]), np.max(positionArray[2, :]) - np.min(positionArray[2, :])]
    # Change lx, ly and lz in the dictionary
    for axis, key in enumerate(systemProps):
        thisMoleculeDict[key] = boxDims[axis]
    # Update the position with the unwrapped_position and then delete that key
    for atomID, position in enumerate(thisMoleculeDict['position']):
        thisMoleculeDict['position'][atomID] = thisMoleculeDict['unwrapped_position'][atomID]
        thisMoleculeDict['image'][atomID] = [0, 0, 0]
    thisMoleculeDict.pop('unwrapped_position')
    # Now centre the molecule based on the average position of its constituent atoms
    thisMoleculeDict = hf.centre(thisMoleculeDict, centre)
    # Finally, update the number of atoms we have in this molecule
    thisMoleculeDict['natoms'] = len(thisMoleculeDict['type'])
    return thisMoleculeDict


def execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList):
    # Main execution loop for the extractMol routine
    moleculeMaster = []
    # Create the moleculeMaster, where each entry is a list of all the AAIDs
    # for the atoms in each molecule
    for molecule in CGToAAIDMaster:
        moleculeAAIDs = []
        for CGSite in list(molecule.keys()):
            moleculeAAIDs += molecule[CGSite][1]
        moleculeMaster.append(sorted(moleculeAAIDs))
    # Iterate through the moleculeMaster, creating new XML files for each
    # molecule in the morphology
    for moleculeNo, moleculeAAIDs in enumerate(moleculeMaster):
        # Use the moleculeAAIDs for each molecule to create a writeable molecule dictionary
        moleculeDict = obtainMoleculeDict(AAMorphologyDict, moleculeAAIDs)
        # Write the molecule dictionary
        hf.writeMorphologyXML(moleculeDict, parameterDict['outputMorphDir'] + '/' + parameterDict['morphology'][:-4] + '/molecules/mol_%04d.xml' % (moleculeNo))


if __name__ == "__main__":
    try:
        pickleFile = sys.argv[1]
    except:
        print("Please specify the pickle file to load to continue the pipeline from this point.")
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = hf.loadPickle(pickleFile)
    execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList)
