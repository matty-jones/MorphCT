import os
import sys
import numpy as np
import cPickle as pickle
import helperFunctions
import csv
import time as T
import random as R
import matplotlib.pyplot as plt


elementaryCharge = 1.60217657E-19 # C
kB = 1.3806488E-23 # m^{2} kg s^{-2} K^{-1}


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
        #### NON-ZERO TRANSFER INTEGRAL ####
        #if pair[6] != 0.0:
        ####################################
        ## LINEARLY SCALED RANDOM NUMBERS ##
        #if R.random() < pair[6]/float(maximumTI):
        ####################################
        #### EXPONENTIAL BOLTZMANN TERM ####
        if R.random() > np.exp(-pair[6]*elementaryCharge/(kB*effectiveT)):
        ####################################
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


def createTCLScript(morphologyName, clusterCommands, highlightClusters):
    xmlFileName = os.getcwd()+'/outputFiles/'+morphologyName+'/morphology/relaxed_'+morphologyName+'.xml'
    hyphenLocs = helperFunctions.findIndex(morphologyName, '-')
    tempName = morphologyName[hyphenLocs[3]+1:hyphenLocs[4]]
    tclLinesToWrite = ['mol delrep 0 0;'] # Load the tcl environment and delete the original representation
    # Now wrap the box and reset the view
    tclLinesToWrite += ['pbc wrap -center origin;', 'pbc box -color black -center origin -width 6;', 'display resetview;']
    # Create the new `faded' material so that we can do cluster highlighting
    tclLinesToWrite += ['material add copy AOEdgy;', 'material rename Material23 Faded;', 'material change opacity Faded 0.02;']
    # The pink looks too similar to the red so change it
    tclLinesToWrite += ['color change rgb 9 1.0 0.29 0.5;']
    # Make all of the atoms faded to begin with (excluding the ones we're going to highlight because otherwise this command dominates in the snapshot)
    if len(highlightClusters) > 0:
        commandsToWrite = [['all and not index']]
    else:
        commandsToWrite = [['all']]
    for cluster in highlightClusters:
        commandsToWrite.append(clusterCommands[cluster])
        commandsToWrite[0] += clusterCommands[cluster][1:]
    for repNo, command in enumerate(commandsToWrite):
        if repNo == 0:
            tclLinesToWrite += ['mol color ColorID '+str(repNo%33)+';', 'mol representation VDW 1.0 8.0;', 'mol material Faded;', 'mol addrep 0;']
            tclLinesToWrite += ['mol modselect '+str(repNo)+' 0 '+' '.join(command)+' and not type H1 C3 C4 C5 C6 C7 C8;']
            continue
        tclLinesToWrite += ['mol color ColorID '+str(repNo%33)+';', 'mol representation VDW 1.0 8.0;', 'mol material AOEdgy;', 'mol addrep 0;']
        tclLinesToWrite += ['mol modselect '+str(repNo)+' 0 '+' '.join(command)+' and not type H1 C3 C4 C5 C6 C7 C8;']


    
    # # Too many reps causes mem error, so have everything faded then highlight the big clusters
    # tclLinesToWrite += ['mol color ColorID 0\n', 'mol representation VDW 1.0 8.0\n', 'mol material AOEdgy\n', 'mol addrep 0\n']
    # for repNo, cluster in enumerate(clusterCommands):
    #     # Create a new representation, coloured by the repNo
    #     # tclLinesToWrite += ['mol color ColorID '+str(repNo%33)+'\n', 'mol representation VDW 1.0 8.0\n', 'mol selection all\n', 'mol material AOEdgy\n', 'mol addrep 0\n']
    #     if repNo in highlightCluster:
    #         tclLinesToWrite += ['mol color ColorID '+str(repNo%33)+'\n', 'mol representation VDW 1.0 8.0\n', 'mol material AOEdgy\n', 'mol addrep 0\n']
    #         tclLinesToWrite += ['mol modselect '+str(repNo)+' 0 '+' '.join(cluster)+' and not type H1 S1 C3 C4 C5 C6 C7 C8\n']

    #     # Too many representations causes memory error
    #     # else:
    #     #     tclLinesToWrite += ['mol color ColorID '+str(repNo%33)+'\n', 'mol representation VDW 1.0 8.0\n', 'mol material Faded\n', 'mol addrep 0\n']
    #     # Now change what we're showing for this repNo and remove the sidechains
    tclFileName = './VMD_'+tempName+'.tcl'
    with open(tclFileName, 'w+') as tclFile:
        tclFile.writelines(tclLinesToWrite)
    print 'TCL file written to', tclFileName
    return float(tempName[1:])



def plotClusterDist(clusterDist, temperature, figSaveDir):
    if len(clusterDist) == 1:
        print "Skipping cluster histogram as only a single cluster"
        return clusterDist[-1], len(clusterDist)
    plt.figure()
    plt.hist(clusterDist)
    plt.xlabel('Size of cluster')
    plt.ylabel('Frequency')
    plt.savefig(figSaveDir+'/ClusterDistT'+str(temperature)+'.png')
    plt.close()
    return clusterDist[-1], len(clusterDist)



def execute(morphologyFile, morphologyName, AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize, chromoDict, singlesData, pairsData):
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
    # print "-= VMD ATOM SELECTIONS PER CLUSTER =-"
    # for command in VMDCommands:
    #     print ' '.join(command)
    # Create a tcl script
    clusterIDs = range(len(atomIDsByCluster.keys()))
    clusterLengths = [len(cluster) for cluster in atomIDsByCluster.values()]
    sortedLengths, sortedIDs = helperFunctions.parallelSort(clusterLengths, clusterIDs)
    # Highlight the second biggest cluster and then fade out all the others
    # This breaks everything because I have too many small clusters.
    # How about I just pick some clusters and highlight those, make the rest into one faded blob.
    # Ten biggest
    highlightClusters = sortedIDs[-10:]
    # All but the biggest
    #highlightClusters = sortedIDs[:-1]
    temperature = createTCLScript(morphologyName, VMDCommands, highlightClusters)

    return sortedLengths, temperature
    # return morphologyFile, AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize, chromoDict, singlesData, pairsData




def loadData(morphologyFile, figSaveDir):
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
    # morphologyFile, AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize, chromoDict, singlesData, pairsData = execute(morphologyFile, morphologyName, AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize, chromoDict, singlesData, pairsData)
    clusterDist, temperature = execute(morphologyFile, morphologyName, AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize, chromoDict, singlesData, pairsData)
    #print "Cluster Distribution =", clusterDist
    sizeOfBiggestCluster, numberOfClusters = plotClusterDist(clusterDist, temperature, figSaveDir)
    return sizeOfBiggestCluster, numberOfClusters, temperature, clusterDist
    # return morphologyFile, AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize, chromoDict, singlesData, pairsData


if __name__ == '__main__':
    sys.setrecursionlimit(10000) # If I use the 01mon system, I hit the recursion limit even though the program is working fine.
    R.seed(32)
    morphologyDir = './outputFiles'

    #for tempVal in np.arange(1000, 10001, 1000):
    #for tempVal in [290, 390, 590, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]:
    tempVal = 3000
    global effectiveT
    effectiveT = tempVal
    tempString = str(tempVal)
    while len(tempString) < 5:
        tempString = '0'+tempString
    saveDir = './percolationGraphs/'+tempString+'K'

    sizeOfBiggestCluster = []
    temperature = []
    numberOfClusters = []
    clusterDists = []
    for morphologyFile in os.listdir(morphologyDir):
        biggestClusterSize, clusterQuantity, temp, clusterDist = loadData(morphologyDir+'/'+morphologyFile, saveDir)
        sizeOfBiggestCluster.append(biggestClusterSize)
        numberOfClusters.append(clusterQuantity)
        temperature.append(temp)
        clusterDists.append(clusterDist)

    plt.figure()
    plt.plot(temperature, sizeOfBiggestCluster)
    plt.xlabel('Temperature')
    plt.ylabel('Size of Biggest Cluster')
    #plt.ylim([0, 1000])
    plt.savefig(saveDir+'/clusterSize.png')
    print "Figure saved to "+saveDir+"/clusterSize.png"

    plt.clf()
    plt.plot(temperature, numberOfClusters)
    plt.xlabel('Temperature')
    plt.ylabel('Number of Clusters')
    plt.savefig(saveDir+'/clusterQuantity.png')
    plt.close()
    print "Figure saved to "+saveDir+"/clusterQuantity.png"

    with open(saveDir+'/clusterData.csv', 'w+') as csvFile:
        csvWriter = csv.writer(csvFile)
        # Data in structure:
        # Temperature, Size of Biggest, Number of Clusters, ClusterDist
        for morphNo in range(len(temperature)):
            csvWriter.writerow([temperature[morphNo], sizeOfBiggestCluster[morphNo], numberOfClusters[morphNo], str(clusterDists[morphNo])])
