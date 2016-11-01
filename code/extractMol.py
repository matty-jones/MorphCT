from hoomd_script import *
import numpy as np
import copy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cPickle as pickle
import helperFunctions
import sys









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
            if sum([x in constraint[1:] for x in moleculeAAIDs]) == len(constraint)-1:
                newConstraint = copy.deepcopy(constraint)
                # Update all of the AAIDs to their new IDs within this molecule
                for location, atomID in enumerate(newConstraint):
                    if location == 0:
                        continue
                    newConstraint[location] = AAIDConversion[str(atomID)]
                # Now add them to the molecule dict
                thisMoleculeDict[constraintType].append(newConstraint)
    # Now need to make sure that everything is unwrapped
    thisMoleculeDict = helperFunctions.addUnwrappedPositions(thisMoleculeDict)
    # Find the (geometric, in case masses aren't specified) centre of the molecule
    positionArray = np.array(thisMoleculeDict['unwrapped_position'])
    centre = [np.average(positionArray[0,:]), np.average(positionArray[1,:]), np.average(positionArray[2,:])]
    # Also find the extent of the molecule so we can change lx, ly and lz
    boxDims = [np.max(positionArray[0,:]) - np.min(positionArray[0,:]), np.max(positionArray[1,:]) - np.min(positionArray[1,:]), np.max(positionArray[2,:]) - np.min(positionArray[2,:])]
    # Change lx, ly and lz in the dictionary
    for axis, key in enumerate(systemProps):
        thisMoleculeDict[key] = boxDims[axis]
    # Update the position with the unwrapped_position and then delete that key
    for atomID, position in enumerate(thisMoleculeDict['position']):
        thisMoleculeDict['position'][atomID] = thisMoleculeDict['unwrapped_position'][atomID]
        thisMoleculeDict['image'][atomID] = [0, 0, 0]
    thisMoleculeDict.pop('unwrapped_position')
    # Now centre the molecule based on the average position of its constituent atoms
    thisMoleculeDict = helperFunctions.centre(thisMoleculeDict, centre)
    # Finally, update the number of atoms we have in this molecule
    thisMoleculeDict['natoms'] = len(thisMoleculeDict['type'])
    return thisMoleculeDict


def execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict):
    moleculeMaster = []
    # Create the moleculeMaster, where each entry is a list of all the AAIDs for the atoms in each molecule
    for molecule in CGToAAIDMaster:
        moleculeAAIDs = []
        for CGSite in molecule.keys():
            moleculeAAIDs += molecule[CGSite][1]
        moleculeMaster.append(sorted(moleculeAAIDs))
    for moleculeNo, moleculeAAIDs in enumerate(moleculeMaster):
        # Use the moleculeAAIDs for each molecule to create a writeable molecule dictionary
        moleculeDict = obtainMoleculeDict(AAMorphologyDict, moleculeAAIDs)
        # Write the molecule dictionary
        helperFunctions.writeMorphologyXML(moleculeDict, parameterDict['outputDir']+'/'+parameterDict['morphology'][:-4]+'/molecules/mol_%04d.xml' % (moleculeNo))
        break







if __name__ == "__main__":
    try:
        pickleFile = sys.argv[1]
    except:
        print "Please specify the pickle file to load to continue the pipeline from this point."
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict = helperFunctions.loadPickle(pickleFile)
    execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict)







#
#
#
#def execute(morphologyFile, AAfileName, inputCGMorphologyDict, inputAAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize):
#    morphologyName = morphologyFile[helperFunctions.findIndex(morphologyFile,'/')[-1]+1:]
#    outputDir = './outputFiles'
#    morphologyList = os.listdir(outputDir)
#    for allMorphologies in morphologyList:
#        if morphologyName in allMorphologies:
#            outputDir += '/'+morphologyName
#            break
#    fileList = os.listdir(outputDir+'/morphology')
#    secondDirList = os.listdir(outputDir)
#    moleculePOSCARS = []
#    if 'molecules' not in secondDirList:
#        os.makedirs(outputDir+'/molecules')
#    else:
#        moleculeFiles = os.listdir(outputDir+'/molecules')
#        for moleculeFile in moleculeFiles:
#            if ('.POSCAR' in moleculeFile) or ('.poscar' in moleculeFile):
#                moleculePOSCARS.append(moleculeFile)
#    moleculeAAIDs = helperFunctions.getAAIDsByMolecule(CGtoAAIDs)
#    if len(moleculePOSCARS) != 0:
#        if len(moleculePOSCARS) == len(moleculeAAIDs):
#            print "All molecule files already treated. Please delete the .POSCAR files to run the set again."
#            return morphologyFile, AAfileName, inputCGMorphologyDict, inputAAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize
#    # GET sSCALE FROM SOMEWHERE ELSE RATHER THAN HARDCODING IT IN HERE!
#    inverseSScale = helperFunctions.getsScale(outputDir, morphologyName)
#    
#    slashList = helperFunctions.findIndex(AAfileName, '/')
#    inputFileName = AAfileName[:slashList[-1]+1]+'relaxed_'+AAfileName[slashList[-1]+1:]
#    print "Loading morphology data..."
#    AAMorphologyDict = helperFunctions.loadMorphologyXML(inputFileName)
#
#    # Unwrap positions
#    AAMorphologyDict = helperFunctions.addUnwrappedPositions(AAMorphologyDict)
#
#    for moleculeNo, AAIDs in enumerate(moleculeAAIDs):
#    # moleculeNo = 63
#    # AAIDs = moleculeAAIDs[63]
#    # Now make the DFT input files
#        atomIDOffset = -np.min(AAIDs)
#        minimumIndex = np.min(AAIDs)
#        maximumIndex = np.max(AAIDs)
#        moleculeDictionary = {'position':[], 'type':[], 'bond':[]}
#        for bond in AAMorphologyDict['bond']:
#            # Can reduce the number of calculations by assuming that the fine-grainer always builds molecules sequentially so the atom IDs go from minimumIndex -> maximumIndex with no breaks
#            if ((bond[1] >= minimumIndex) and (bond[1] <= maximumIndex)) or ((bond[2] >= minimumIndex) and (bond[2] <= maximumIndex)):
#                moleculeDictionary['bond'].append([bond[0], bond[1]+atomIDOffset, bond[2]+atomIDOffset])
#        molID = str(moleculeNo)
#        while len(molID) < 3:
#            molID = '0'+molID
#        poscarName = outputDir+'/molecules/mol'+molID+'.POSCAR'
#        nAtoms = 0
#        for AAID in AAIDs:
#            moleculeDictionary['position'].append(AAMorphologyDict['unwrapped_position'][AAID])
#            moleculeDictionary['type'].append(AAMorphologyDict['type'][AAID])
#            nAtoms += 1
#        for key in ['lx', 'ly', 'lz']:
#            moleculeDictionary[key] = AAMorphologyDict[key]
#        moleculeDictionary['natoms'] = nAtoms
#        moleculeDictionary = helperFunctions.addMasses(moleculeDictionary)
#        moleculeCOM = helperFunctions.calcCOM(moleculeDictionary['position'], moleculeDictionary['mass'])
#        moleculeDictionary = helperFunctions.centre(moleculeDictionary, moleculeCOM)
#        moleculeDictionary = helperFunctions.scale(moleculeDictionary, inverseSScale)
#        moleculeDictionary = helperFunctions.alignMolecule(moleculeDictionary, [0,1,0]) # Align along y-axis
#        helperFunctions.writePOSCARFile(moleculeDictionary, poscarName)
#        # #Now need to write the pickle again now that we have calculated moleculeAAIDs
#        # pickleFileName = './outputFiles/'+morphologyName+'/morphology/'+morphologyName+'.pickle'
#        # print "Updating pickle file..."
#        # with open(pickleFileName, 'w+') as pickleFile:
#        #     pickle.dump((AAfileName, inputCGMorphologyDict, inputAAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize), pickleFile)
#        # print "Pickle file written to", pickleFileName
#    return morphologyFile, AAfileName, inputCGMorphologyDict, inputAAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize
#    
#
#def loadPickle(morphologyFile):
#    morphologyName = morphologyFile[helperFunctions.findIndex(morphologyFile,'/')[-1]+1:]
#    outputDir = './outputFiles'
#    morphologyList = os.listdir(outputDir)
#    for allMorphologies in morphologyList:
#        if morphologyName in allMorphologies:
#            outputDir += '/'+morphologyName
#            break
#    fileList = os.listdir(outputDir+'/morphology')
#    pickleFound = False
#    for fileName in fileList:
#        if fileName == morphologyName+'.pickle':
#            pickleLoc = outputDir+'/morphology/'+fileName
#            pickleFound = True
#    if pickleFound == False:
#        print "Pickle file not found. Please run morphCT.py again to create the required HOOMD inputs."
#        exit()
#    print "Pickle found at", str(pickleLoc)+"."
#    print "Loading atom data..."
#    with open(pickleLoc, 'r') as pickleFile:
#        (AAfileName, inputCGMorphologyDict, inputAAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize) = pickle.load(pickleFile)
#    morphologyFile, AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize = execute(morphologyFile, AAfileName, inputCGMorphologyDict, inputAAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize)
#    return morphologyFile, AAfileName, inputCGMorphologyDict, inputAAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize
#
#
#if __name__ == "__main__":
#    morphologyFile = sys.argv[1]
#    loadPickle(morphologyFile)
