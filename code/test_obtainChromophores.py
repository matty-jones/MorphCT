import pytest
import obtainChromophores
import executeZINDO
import helperFunctions
import copy
import random as R

def testFindNeighbours(pickleFile):
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, oldChromophoreList = helperFunctions.loadPickle(pickleFile)
    emptyCutOffChromophoreList = copy.deepcopy(oldChromophoreList)
    emptyVoronoiChromophoreList = copy.deepcopy(oldChromophoreList)
    for chromoIndex, chromo in enumerate(oldChromophoreList):
        emptyCutOffChromophoreList[chromoIndex].neighbours = []
        emptyCutOffChromophoreList[chromoIndex].dissociationNeighbours = []
        emptyVoronoiChromophoreList[chromoIndex].neighbours = []
        emptyVoronoiChromophoreList[chromoIndex].dissociationNeighbours = []
    simDims = [[-axis/2.0, axis/2.0] for axis in [AAMorphologyDict[boxLength] for boxLength in ['lx', 'ly', 'lz']]]
    parameterDict['maximumHoleHopDistance'] = 14.2
    parameterDict['maximumElectronHopDistance'] = 14.2
    oldChromophoreList = obtainChromophores.determineNeighboursCutOff(emptyCutOffChromophoreList, parameterDict, simDims)
    newChromophoreList = obtainChromophores.determineNeighboursVoronoi(emptyVoronoiChromophoreList, parameterDict, simDims)
    for listName in [oldChromophoreList, newChromophoreList]:
        chromoID = 653
        print(chromoID)
        print(' '.join(list(map(str, listName[chromoID].AAIDs + [item for sublist in [listName[x[0]].AAIDs for x in listName[chromoID].neighbours] for item in sublist]))))
        print(' '.join(list(map(str, listName[chromoID].AAIDs + [item for sublist in [listName[x[0]].AAIDs for x in listName[chromoID].dissociationNeighbours] for item in sublist]))) + '\n')


def testWriteORCAOutput(pickleFile):
    # One of the chromophores in the corner is #1198
    R.seed(8585)
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = helperFunctions.loadPickle(pickleFile)
    for chromo in chromophoreList:
        chromo.neighbours = []
        chromo.dissociationNeighbours = []
    simDims = [[-axis/2.0, axis/2.0] for axis in [AAMorphologyDict[boxLength] for boxLength in ['lx', 'ly', 'lz']]]
    chromophoreList = obtainChromophores.determineNeighboursVoronoi(chromophoreList, parameterDict, simDims)
    #for runNumber in range(20):
    runNumber = 0
    chromoID = 1187#R.randint(0, len(chromophoreList))
    chromophore1 = chromophoreList[chromoID]
    print("\n\n", chromoID, chromophore1.neighbours)
    AAIDs = chromophore1.AAIDs
    images = [[0, 0, 0] for i in range(len(chromophore1.AAIDs))]
    for index, neighbourChromo in enumerate(chromophore1.neighbours):
        chromophore2 = chromophoreList[neighbourChromo[0]]
        AAIDs += chromophore2.AAIDs
        images += [neighbourChromo[1] for i in range(len(chromophore2.AAIDs))]
    executeZINDO.writeOrcaInp(AAMorphologyDict, AAIDs, images, None, None, './testAssets/outputFiles/testORCAInput_%03d.inp' % (runNumber))



if __name__ == "__main__":
    pickleFile = 'testAssets/bilayerBCC/code/bilayerBCC.pickle'
    testFindNeighbours(pickleFile)
    #testWriteORCAOutput(pickleFile)
