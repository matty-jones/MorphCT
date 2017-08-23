import sys
import os
import helperFunctions
import cPickle as pickle
import copy


def loadChromoPickle(path):
    print "Loading", str(path)+"..."
    with open(path, 'r') as pickleFile:
        chromoDict = pickle.load(pickleFile)
    return chromoDict


def loadMorphPickle(path):
    print "Loading", str(path)+"..."
    with open(path, 'r') as pickleFile:
        AAfileName, CGMoleculeDict, UnrelaxedAAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize = pickle.load(pickleFile)
    return AAfileName, CGMoleculeDict, UnrelaxedAAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize


def getChromoID(fileName):
    slashList = helperFunctions.findIndex(fileName, '/')
    if slashList != None:
        fileName = fileName[slashList[-1]+1:]
    return map(int, ('_'+fileName[:-4]).split('_chromo')[1:])


def getChromosToKeep(dirName):
    fileNames = os.listdir(dirName)
    chromosToKeep = []
    for fileName in fileNames:
        chromoIDs = getChromoID(fileName)
        for chromoID in chromoIDs:
            if int(chromoID) not in chromosToKeep:
                chromosToKeep.append(int(chromoID))
    return chromosToKeep

def getAtomIDsToKeep(chromosToKeep, chromoDict):
    moleculesToKeep = []
    for chromophore in chromoDict.keys():
        if chromoDict[chromophore]['periodic'] == False:
            if chromophore in chromosToKeep:
                moleculesToKeep.append(chromoDict[chromophore]['molID'])
    moleculesToKeep = list(set(moleculesToKeep))
    atomIDsToKeep = []
    for chromoID in chromoDict.keys():
        if chromoDict[chromoID]['periodic'] == False:
            if chromoDict[chromoID]['molID'] in moleculesToKeep:
                atomIDsToKeep += chromoDict[chromoID]['atomID']
    return sorted(atomIDsToKeep), sorted(moleculesToKeep)


def trimCGtoAAIDs(CGtoAAIDs, moleculeIDs, AAIDLookup):
    newCGtoAAIDs = []
    for moleculeID in moleculeIDs:
        for key in CGtoAAIDs[moleculeID].keys():
            newAAIDs = []
            for oldAAID in CGtoAAIDs[moleculeID][key][1]:
                newAAIDs.append(AAIDLookup[oldAAID])
            CGtoAAIDs[moleculeID][key][1] = newAAIDs
        newCGtoAAIDs.append(CGtoAAIDs[moleculeID])
    return newCGtoAAIDs


def trimCGMoleculeDict(CGMoleculeDict, CGtoAAIDs):
    CGIDs = []
    for molecule in CGtoAAIDs:
        CGIDs += molecule.keys()
    lookup = {}
    for newCGID, oldCGID in enumerate(CGIDs):
        lookup[oldCGID] = newCGID
    newCGMoleculeDict = copy.deepcopy(CGMoleculeDict)
    newCGMoleculeDict['natoms'] = len(CGIDs)
    for key in ['body', 'diameter', 'image', 'position', 'unwrapped_position', 'mass', 'velocity', 'type']:
            print "Updating newCGMoleculeDict["+key+"]..."
            newData = []
            for CGID in CGIDs:
                newData.append(CGMoleculeDict[key][CGID])
            newCGMoleculeDict[key] = newData
    for key in ['bond']:
            print "Updating newCGMoleculeDict["+key+"] (may take a while)..."
            newData = []
            for element in CGMoleculeDict[key]:
                important = False
                for CGID in element[1:]:
                    if CGID in CGIDs:
                        important = True
                        break
                    # else:
                    #     important = False
                    #     break
                if important == True:
                    for index, value in enumerate(element):
                        if index == 0:
                            continue
                        #element[index] = lookup[value]
                    newData.append(element)
            newCGMoleculeDict[key] = newData
    return CGMoleculeDict


def trimBonds(AAMorphologyDict, AAIDLookup):
    print len(AAMorphologyDict['bond'])
    print "Fixing AAMorphologyDict bonds..."
    newMorphologyBonds = []
    for bondno, bond in enumerate(AAMorphologyDict['bond']):
        if (bond[1] in AAIDLookup) and (bond[2] in AAIDLookup):
            newMorphologyBonds.append([bond[0], AAIDLookup[bond[1]], AAIDLookup[bond[2]]])
    AAMorphologyDict['bond'] = newMorphologyBonds
    print len(AAMorphologyDict['bond'])
    return AAMorphologyDict


if __name__ == "__main__":
    generateNewMorphologyXML = False
    generateNewPickleFile = True

    chromoDict = loadChromoPickle('./chromophores.pickle')
    morphologyDict = helperFunctions.loadMorphologyXML('./relaxed_p1-L15-f0.0-P0.1-T2.25-e0.5.xml')#, sigma=3.0)
    chromosToKeep = getChromosToKeep('./brokenInps')
    atomIDsToKeep, moleculesToKeep = getAtomIDsToKeep(chromosToKeep, chromoDict)


    lookup = {}
    for newAAID, oldAAID in enumerate(atomIDsToKeep):
        lookup[oldAAID] = newAAID


    if generateNewMorphologyXML == True:
        newMorphDict = copy.deepcopy(morphologyDict)
        newMorphDict['natoms'] = len(atomIDsToKeep)
        # print lookup
        # print sorted(lookup.values())
        for key in ['body', 'diameter', 'image', 'position', 'charge', 'mass', 'velocity', 'type']:
            print "Updating newMorphDict["+key+"]..."
            newData = []
            for atomID in atomIDsToKeep:
                newData.append(morphologyDict[key][atomID])
            newMorphDict[key] = newData
        for key in ['bond', 'angle', 'dihedral', 'improper']:
            print "Updating newMorphDict["+key+"] (may take a while)..."
            newData = []
            for element in morphologyDict[key]:
                important = False
                for atomID in element[1:]:
                    if atomID in atomIDsToKeep:
                        important = True
                        break
                    # else:
                    #     important = False
                    #     break
                if important == True:
                    for index, value in enumerate(element):
                        if index == 0:
                            continue
                        element[index] = lookup[value]
                    newData.append(element)
            newMorphDict[key] = newData
        helperFunctions.writeMorphologyXML(newMorphDict, './outputFiles/testMorph.xml')

    if generateNewPickleFile == True:
        AAfileName, CGMoleculeDict, UnrelaxedAAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize = loadMorphPickle('./p1-L15-f0.0-P0.1-T2.25-e0.5.pickle')
        AAfileName = './outputFiles/testMorph/morphology/testMorph.xml'
        CGtoAAIDs = trimCGtoAAIDs(CGtoAAIDs, moleculesToKeep, lookup)
        CGMoleculeDict = trimCGMoleculeDict(CGMoleculeDict, CGtoAAIDs)
        UnrelaxedAAMorphologyDict = trimBonds(UnrelaxedAAMorphologyDict, lookup)
        toPickle = (AAfileName, CGMoleculeDict, UnrelaxedAAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize)
        with open('./outputFiles/testMorph.pickle', 'w+') as pickleFile:
            pickle.dump(toPickle, pickleFile)


    # for analyseMolecules need morphologyFile, AAfileName, CGMoleculeDict, CGtoAAIDs, boxSize
    # Don't need UnrelaxedAAMorphologyDict, moleculeAAIDs
