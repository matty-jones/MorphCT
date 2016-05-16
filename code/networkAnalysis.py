import os
import sys
import numpy as np
import cPickle as pickle
import helperFunctions
import csv
import time as T
import random as R


def loadCSVs(outputDir):
    CSVDir = outputDir+'/chromophores'
    singlesData = {}
    pairsData = []
    with open(CSVDir+'/singles.csv', 'r') as singlesFile:
        singlesReader = csv.reader(singlesFile, delimiter=',')
        for row in singlesReader:
            singlesData[int(float(row[0]))] = map(float, row[1:])
    with open(CSVDir+'/pairs.csv', 'r') as pairsFile:
        pairsReader = csv.reader(pairsFile, delimiter=',')
        for row in pairsReader:
            pairsData.append([int(float(row[0])), int(float(row[1]))] + [x for x in map(float, row[2:])])
    return singlesData, pairsData


def getConnectedChromos(pairsData):
    pairsDataArray = np.array(pairsData)
    maximumTI = np.amax(np.array(pairsData)[:,6])
    connectedChromoDict = {}
    for pair in pairsData:
        if pair[0] not in connectedChromoDict:
            connectedChromoDict[pair[0]] = []
        if pair[1] not in connectedChromoDict:
            connectedChromoDict[pair[1]] = []
        # If non-zero transfer integral:
        if pair[6] != 0.0:
            # Add `neighbours' to dict if a random number is < scaled TI
            if R.random() < pair[6]/float(maximumTI):
                connectedChromoDict[pair[0]].append(pair[1])
                connectedChromoDict[pair[1]].append(pair[0])
    for index, neighbours in connectedChromoDict.iteritems():
        connectedChromoDict[index] = sorted(neighbours)
    return connectedChromoDict


def updateClusterList(chromoID, completeNeighbourList, clusterList):
    neighbourList = completeNeighbourList[chromoID]
    for neighbour in neighbourList:
        if clusterList[neighbour] > clusterList[chromoID]:
            previousCluster = clusterList[neighbour]
            clusterList[neighbour] = clusterList[chromoID]
            #print "Current ChromoID =", chromoID, "moving", neighbour, "from", previousCluster, "to", clusterList[neighbour]
            updateClusterList(neighbour, completeNeighbourList, clusterList)
        elif clusterList[neighbour] < clusterList[chromoID]:
            previousCluster = clusterList[neighbour]
            clusterList[chromoID] = clusterList[neighbour]
            #print "Current ChromoID =", chromoID, "moving", chromoID, "from", previousCluster, "to", clusterList[chromoID]
            updateClusterList(chromoID, completeNeighbourList, clusterList)
    return clusterList


def getClusterDict(connectedChromoDict):
    clusterDict = {}
    for chromoID in connectedChromoDict.keys():
        clusterDict[chromoID] = chromoID
    for chromoID in sorted(connectedChromoDict.keys()):
        clusterDict = updateClusterList(chromoID, connectedChromoDict, clusterDict)
    # Flip the dictionary around as we did before
    clustDict = {}
    for chromoID in clusterDict.keys():
        clusterID = clusterDict[chromoID]
        if clusterID not in clustDict:
            clustDict[clusterID] = [chromoID]
        else:
            clustDict[clusterID].append(chromoID)
    return clustDict


def getAtomIDs(clusterDict, chromoDict):
    #print chromoDict[6541]
    atomIDsByCluster = {}
    for clusterID in clusterDict.keys():
        AAIDs = []
        for chromoID in clusterDict[clusterID]:
            AAIDs += chromoDict[chromoID]['atomID']
        atomIDsByCluster[clusterID] = AAIDs
    return atomIDsByCluster


def generateVMDSelection(AAIDList, printFlag):
    # Work out a compact way of VMD taking AAIDList and importing it to tcl
    selectionCommand = ['index']
    AAIDList = sorted(AAIDList)
    previousAtom = None
    inRange = False
    for AAID in AAIDList:
        if AAID - 1 != previousAtom:
            if previousAtom != None:
                selectionCommand.append('to')
                selectionCommand.append(str(previousAtom))
            selectionCommand.append(str(AAID))
        previousAtom = AAID
    if selectionCommand[-1] != str(previousAtom):
        selectionCommand.append('to')
        selectionCommand.append(str(previousAtom))
    if printFlag == True:
        print AAIDList
        print selectionCommand
    return selectionCommand



def execute(morphologyFile, AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize, chromoDict, singlesData, pairsData):
    realChromoIDs = []
    for chromoID in chromoDict.keys():
        if chromoDict[chromoID]['realChromoID'] not in realChromoIDs:
            realChromoIDs.append(chromoDict[chromoID]['realChromoID'])
    connectedChromoDict = getConnectedChromos(pairsData)
    clusterDict = getClusterDict(connectedChromoDict)

    # for chromoID, neighbours in connectedChromoDict.iteritems():
    #     print chromoID, neighbours
    # print clusterDict
    atomIDsByCluster = getAtomIDs(clusterDict, chromoDict)
    # for clusterID, atomsIDs in enumerate(atomIDsByCluster):
    #     print clusterID, atomsIDs
    # print len(clusterDict)
    # print len(atomIDsByCluster)
    VMDCommands = []
    printFlag = False
    for clusterID in atomIDsByCluster.keys():
        selectionCommand = generateVMDSelection(atomIDsByCluster[clusterID], printFlag)
        VMDCommands.append(selectionCommand)
        #printFlag = True
    print "-= VMD ATOM SELECTIONS PER CLUSTER =-"
    for command in VMDCommands:
        print ' '.join(command)
    return morphologyFile, AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize, chromoDict, singlesData, pairsData




def loadData(morphologyFile):
    morphologyName = morphologyFile[helperFunctions.findIndex(morphologyFile,'/')[-1]+1:]
    outputDir = './outputFiles'
    morphologyList = os.listdir(outputDir)
    for allMorphologies in morphologyList:
        if morphologyName in allMorphologies:
            outputDir += '/'+morphologyName
            break
    morphPickleFound = False
    chromoPickleFound = False
    for fileName in os.listdir(outputDir+'/morphology'):
        if fileName == morphologyName+'.pickle':
            morphPickleLoc = outputDir+'/morphology/'+fileName
            morphPickleFound = True
        elif fileName == 'chromophores.pickle':
            chromoPickleLoc = outputDir+'/morphology/'+fileName
            chromoPickleFound = True
    if morphPickleFound == False:
        print "Morphology pickle file not found. Please run morphCT.py again to create the required HOOMD inputs."
        exit()
    if chromoPickleFound == False:
        print "Chromophore pickle file not found. Please run analyseMolecules.py again to create the required chromophore data."
    print "Morphology pickle found at", str(morphPickleLoc)+"."
    print "Loading data..."
    with open(morphPickleLoc, 'r') as pickleFile:
        (AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize) = pickle.load(pickleFile)
    print "Chromophore pickle found at", str(chromoPickleLoc)+"."
    print "Loading data..."
    with open(chromoPickleLoc, 'r') as pickleFile:
        chromoDict = pickle.load(pickleFile)
    singlesData, pairsData = loadCSVs(outputDir)
    morphologyFile, AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize, chromoDict, singlesData, pairsData = execute(morphologyFile, AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize, chromoDict, singlesData, pairsData)
    return morphologyFile, AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize, chromoDict, singlesData, pairsData


if __name__ == '__main__':
    sys.setrecursionlimit(10000) # If I use the 01mon system, I hit the recursion limit even though the program is working fine.
    R.seed(32)
    morphologyFile = sys.argv[1]
    loadData(morphologyFile)
