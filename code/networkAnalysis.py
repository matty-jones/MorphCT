import os
import sys
import numpy as np
import cPickle as pickle
import helperFunctions
import csv
import time as T


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
    connectedChromoDict = {}
    for pair in pairsData:
        if pair[0] not in connectedChromoDict:
            connectedChromoDict[pair[0]] = []
        if pair[1] not in connectedChromoDict:
            connectedChromoDict[pair[1]] = []
        # If non-zero transfer integral:
        if pair[-1] != 0.0:
            # Add `neighbours' to dict
            connectedChromoDict[pair[0]].append(pair[1])
            connectedChromoDict[pair[1]].append(pair[0])
    for index, neighbours in connectedChromoDict.iteritems():
        connectedChromoDict[index] = sorted(neighbours)
    return connectedChromoDict


def getConnectedChromos2(pairsData):
    connectedChromoDict = {}
    for pair in pairsData:
        if pair[0] not in connectedChromoDict:
            connectedChromoDict[pair[0]] = []
        if pair[1] not in connectedChromoDict:
            connectedChromoDict[pair[1]] = []
        # If non-zero transfer integral:
        if pair[-1] != 0.0:
            # Add `neighbours' to dict
            if pair[0] > pair[1]:
                connectedChromoDict[pair[1]].append(pair[0])
            else:
                connectedChromoDict[pair[0]].append(pair[1])
    for index, neighbours in connectedChromoDict.iteritems():
        connectedChromoDict[index] = sorted(neighbours)
    return connectedChromoDict


def getClusterDict(connectedChromoDict):
    clusterDict = {}
    for chromoID in connectedChromoDict.keys():
        clusterDict[chromoID] = chromoID
    while True:
        # Run this chunk of code until it stops updating clusters
        numberOfChanges = 0
        for chromoID, neighbours in connectedChromoDict.iteritems():
            for neighbour in neighbours:
                if clusterDict[chromoID] < clusterDict[neighbour]:
                    if clusterDict[neighbour] != clusterDict[chromoID]:
                        numberOfChanges += 1
                    clusterDict[neighbour] = clusterDict[chromoID]
                else:
                    if clusterDict[chromoID] != clusterDict[neighbour]:
                        numberOfChanges += 1
                    clusterDict[chromoID] = clusterDict[neighbour]
        if numberOfChanges == 0:
            break
    # ClusterDict is currently {ChromoID: Cluster No}.
    # It would be more useful for plotting etc. to flip it around
    clustDict = {}
    for chromoID in clusterDict.keys():
        clusterID = clusterDict[chromoID]
        if clusterID not in clustDict:
            clustDict[clusterID] = [chromoID]
        else:
            clustDict[clusterID].append(chromoID)
    return clustDict


def updateClusterList(chromoID, completeNeighbourList, clusterList):
    neighbourList = completeNeighbourList[chromoID]
    for neighbour in neighbourList:
        if clusterList[neighbour] > clusterList[chromoID]:
            previousCluster = clusterList[neighbour]
            clusterList[neighbour] = clusterList[chromoID]
            if (chromoID == 2) or (chromoID == 113) or (chromoID == 308) or (neighbour == 2) or (neighbour == 113) or (neighbour == 308):
                print "Current ChromoID =", chromoID, "moving", neighbour, "from", previousCluster, "to", clusterList[neighbour]
            updateClusterList(neighbour, completeNeighbourList, clusterList)
        elif clusterList[neighbour] < clusterList[chromoID]:
            previousCluster = clusterList[neighbour]
            clusterList[chromoID] = clusterList[neighbour]
            if (chromoID == 2) or (chromoID == 113) or (chromoID == 308) or (neighbour == 2) or (neighbour == 113) or (neighbour == 308):
                print "Current ChromoID =", chromoID, "moving", chromoID, "from", previousCluster, "to", clusterList[chromoID]
            updateClusterList(chromoID, completeNeighbourList, clusterList)
    return clusterList


def getClusterDict2(connectedChromoDict):
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
    for key in chromoDict.keys():
        print chromoDict[key]['chromoID'], chromoDict[key]['realChromoID']
    atomIDsByCluster = {}
    for clusterID in clusterDict.keys():
        chromoIDs = clusterDict[clusterID]
        AAIDs = []
        for chromoID in chromoIDs:
            AAIDs += chromoDict[chromoID]['atomID']
        atomIDsByCluster[clusterID] = AAIDs
    return atomIDsByCluster


def generateVMDSelection(AAIDList):
    AAIDList = sorted(AAIDList)
    print AAIDList
    exit()



def execute(morphologyFile, AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize, chromoDict, singlesData, pairsData):
    realChromoIDs = []
    for chromoID in chromoDict.keys():
        if chromoDict[chromoID]['realChromoID'] not in realChromoIDs:
            realChromoIDs.append(chromoDict[chromoID]['realChromoID'])
    t0 = T.time()
    connectedChromoDict = getConnectedChromos(pairsData)
    clusterDict = getClusterDict(connectedChromoDict)
    t1 = T.time()
    print "While loop time taken =", t1-t0, "s."
    t2 = T.time()
    connectedChromoDict2 = getConnectedChromos2(pairsData)
    clusterDict2 = getClusterDict2(connectedChromoDict)
    t3 = T.time()
    print "Recursive Fn time taken =", t3-t2, "s."
    
    print "\n\n"
    for chromo, neighbours in connectedChromoDict.iteritems():
        print chromo, neighbours, connectedChromoDict2[chromo]
    print "\n\n"
    
    print clusterDict
    print clusterDict2


    print "Number of clusters in clusterDict =", len(clusterDict.keys())
    print "Number of clusters in clusterDict2 =", len(clusterDict2.keys())
    # check1 = []
    # check2 = []
    # for i, chromos in clusterDict.iteritems():
    #     check1 += chromos
    #     check1.append(i)
    # for i, chromos in clusterDict.iteritems():
    #     check2 += chromos
    #     check2.append(i)
    # print len(check1), len(list(set(check1)))
    # print len(check2), len(list(set(check2)))
    
    if clusterDict2 != clusterDict:
        raise SystemError("Dictionaries not the same.")
    print "Dictionaries identical."
    exit()
    atomIDsByCluster = getAtomIDs(clusterDict, chromoDict)
    VMDCommands = []
    for clusterID in atomIDsByCluster.keys():
        VMDCommand = ['index']
        selectionCommand = generateVMDSelection(atomIDsByCluster[clusterID])
        VMDCommand += selectionCommand
        print VMDCommand
        exit()
        VMDCommands.append(VMDCommand)
        
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
    morphologyFile = sys.argv[1]
    loadData(morphologyFile)
