import pytest
from morphct.code import obtain_chromophores
from morphct.code import execute_ZINDO
from morphct.code import helper_functions as hf
import copy
import random as R
import numpy as np

def testFindNeighbours(pickleFile):
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, oldChromophoreList = hf.loadPickle(pickleFile)
    emptyCutOffChromophoreList = copy.deepcopy(oldChromophoreList)
    emptyVoronoiChromophoreList = copy.deepcopy(oldChromophoreList)
    for chromoIndex, chromo in enumerate(oldChromophoreList):
        emptyCutOffChromophoreList[chromoIndex].neighbours = []
        emptyCutOffChromophoreList[chromoIndex].dissociationNeighbours = []
        emptyVoronoiChromophoreList[chromoIndex].neighbours = []
        emptyVoronoiChromophoreList[chromoIndex].dissociationNeighbours = []
    simDims = [[-axis/2.0, axis/2.0] for axis in [AAMorphologyDict[boxLength] for boxLength in ['lx', 'ly', 'lz']]]
    parameterDict['maximumHoleHopDistance'] = 10.0
    parameterDict['maximumElectronHopDistance'] = 10.0
    oldChromophoreList = obtain_chromophores.determineNeighboursCutOff(emptyCutOffChromophoreList, parameterDict, simDims)
    newChromophoreList = obtain_chromophores.determineNeighboursVoronoi(emptyVoronoiChromophoreList, parameterDict, simDims)
    for listName in [oldChromophoreList, newChromophoreList]:
        chromoID = 653
        print(chromoID)
        print(' '.join(list(map(str, listName[chromoID].AAIDs + [item for sublist in [listName[x[0]].AAIDs for x in listName[chromoID].neighbours] for item in sublist]))))
        print(' '.join(list(map(str, listName[chromoID].AAIDs + [item for sublist in [listName[x[0]].AAIDs for x in listName[chromoID].dissociationNeighbours] for item in sublist]))) + '\n')


def testWriteORCAOutput(pickleFile):
    # One of the chromophores in the corner is #1198
    R.seed(8585)
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = hf.loadPickle(pickleFile)
    for chromo in chromophoreList:
        chromo.neighbours = []
        chromo.dissociationNeighbours = []
    simDims = [[-axis/2.0, axis/2.0] for axis in [AAMorphologyDict[boxLength] for boxLength in ['lx', 'ly', 'lz']]]
    chromophoreList = obtain_chromophores.determineNeighboursCutOff(chromophoreList, parameterDict, simDims)
    #chromophoreList = obtain_chromophores.determineNeighboursVoronoi(chromophoreList, parameterDict, simDims)
    #for runNumber in range(20):
    parameterDict['outputMorphDir'] = './testAssets/outputFiles'
    parameterDict['morphology'] = ''
    execute_ZINDO.createInputFiles(chromophoreList, AAMorphologyDict, parameterDict)
    #chromoID = 2487#R.randint(0, len(chromophoreList))
    #chromophore1 = chromophoreList[chromoID]
    #print("\n\n", chromoID, chromophore1.neighbours)
    #for index, neighbourChromo in enumerate(chromophore1.neighbours):
    #    AAIDs = chromophore1.AAIDs
    #    images = [[0, 0, 0] for i in range(len(chromophore1.AAIDs))]
    #    chromophore2 = chromophoreList[neighbourChromo[0]]
    #    AAIDs += chromophore2.AAIDs
    #    images += [neighbourChromo[1] for i in range(len(chromophore2.AAIDs))]
    #    execute_ZINDO.writeOrcaInp(AAMorphologyDict, AAIDs, images, None, None, './testAssets/outputFiles/testORCAInput_%03d.inp' % (index))


def testCheckPeriodicNeighbours(pickleFile):
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = hf.loadPickle(pickleFile)
    for chromo in chromophoreList:
        chromo.neighbours = []
        chromo.dissociationNeighbours = []
    simDims = [[-axis/2.0, axis/2.0] for axis in [AAMorphologyDict[boxLength] for boxLength in ['lx', 'ly', 'lz']]]
    chromophoreList = obtain_chromophores.determineNeighboursVoronoi(chromophoreList, parameterDict, simDims)
    chromoID = R.randint(0, len(chromophoreList))
    print(chromophoreList[chromoID].neighbours)
    print("\nOriginal =", ' '.join(map(str, chromophoreList[chromoID].AAIDs)))
    neighbour1String = "In-image neighbours = "
    neighbour2String = "Out-of-image neighbours = "
    for [neighbourID, image] in chromophoreList[chromoID].neighbours:
        if np.array_equal(image, [0, 0, 0]):
            neighbour1String += ' '.join(map(str, chromophoreList[neighbourID].AAIDs)) + ' '
        else:
            neighbour2String += ' '.join(map(str, chromophoreList[neighbourID].AAIDs)) + ' '
    print("\n")
    print(neighbour1String)
    print("\n")
    print(neighbour2String)




if __name__ == "__main__":
    #pickleFile = 'testAssets/bilayerBCC/code/bilayerBCC.pickle'
    pickleFile = 'testAssets/p3ht/code/p1-L15-f0.0-P0.1-T1.5-e0.5.pickle'
    #testFindNeighbours(pickleFile)
    testWriteORCAOutput(pickleFile)
    #testPeriodicNeighbours(pickleFile)
