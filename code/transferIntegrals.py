import numpy as np
import sys
import os
import helperFunctions
import csv
import subprocess as sp
import multiprocessing as mp
import cPickle as pickle
import time as T


class ORCAError(Exception):
    def __init__(self, fileName):
        self.string = "No molecular orbital data present for "+str(fileName)
    def __str__(self):
        return self.string
        
class chromophore:
    def __init__(self, inputFile, chromoID, singleChromos = None, chromoData = None):
        if chromoData != None:
            if len(chromoData) == 9:
                # Data for Single
                self.chromo1ID = int(float(chromoData[0]))
                self.position = np.array([float(chromoData[1]), float(chromoData[2]), float(chromoData[3])])
                self.HOMO_1 = float(chromoData[4])
                self.HOMO = float(chromoData[5])
                self.LUMO = float(chromoData[6])
                self.LUMO_1 = float(chromoData[7])
                self.chromoLength = int(float(chromoData[8]))
                self.error = 0
                return
            elif len(chromoData) == 7:
                # Data for Pair
                self.chromo1ID = int(float(chromoData[0]))
                self.chromo2ID = int(float(chromoData[1]))
                self.HOMO_1 = float(chromoData[2])
                self.HOMO = float(chromoData[3])
                self.LUMO = float(chromoData[4])
                self.LUMO_1 = float(chromoData[5])
                self.Tij = float(chromoData[6])
                self.error = 0
                return
            else:
                raise SystemError('Unexpected read chromoData Length')
        self.inputFile = inputFile
        self.chromo1ID = int(chromoID[0])
        if len(chromoID) == 1:
            self.ORCAType = 'single'
        else:
            self.ORCAType = 'pair'            
        try:
            self.loadORCAOutput(inputFile)
        except ORCAError as errorMessage:
            print errorMessage
            self.error = 1
            return
        if len(chromoID) == 2:
            self.chromo2ID = int(chromoID[1])
            self.Tij = self.calculateTransferIntegral(singleChromos)
        self.error = 0


    def loadORCAOutput(self, fileName):
        with open(fileName, 'r') as orcaFile:
            dataFile = orcaFile.readlines()
        recordPosData = False
        recordMOData = False
        listOfPositions = []
        listOfMasses = []
        orbitalData = []
        self.chromoLength = 0
        for line in dataFile:
            if 'CARTESIAN COORDINATES (ANGSTROEM)' in line:
                recordPosData = True
                continue
            if recordPosData == True:
                if 'CARTESIAN COORDINATES (A.U.)' in line:
                    recordPosData = False
                for element in line[:-1].split(' '):
                    if len(element) != 0:
                        try:
                            coordinate = float(element)
                            listOfPositions[-1].append(coordinate)
                        except:
                            if len(element) <= 2:
                                listOfPositions.append([])
                                if element == 'H':
                                    listOfMasses.append(1.00794)
                                elif element == 'C':
                                    listOfMasses.append(12.0107)
                                elif element == 'S':
                                    if self.ORCAType == 'single':
                                        self.chromoLength += 1
                                    listOfMasses.append(32.0660)
                                else:
                                    print element
                                    raise SystemError('Unknown element')
            if 'ORBITAL ENERGIES' in line:
                recordMOData = True
                continue
            if recordMOData == True:
                if 'MOLECULAR ORBITALS' in line:
                    # Don't need anything else from the output file
                    break
                dataInLine = []
                for element in line.split(' '):
                    try:
                        dataInLine.append(float(element))
                    except:
                        pass
                if len(dataInLine) == 4:
                    orbitalData.append(dataInLine)
        self.position = helperFunctions.calcCOM(listOfPositions, listOfMasses)
        for i in range(len(orbitalData)):
            if orbitalData[i][1] == 0:
                # This line is the first unoccupied orbital - i.e. LUMO
                self.LUMO = orbitalData[i][3]
                self.HOMO = orbitalData[i-1][3]
                self.HOMO_1 = orbitalData[i-2][3]
                self.LUMO_1 = orbitalData[i+1][3]
                # Don't need any other orbitals
                break
        if recordMOData == False:
            # Molecular orbital data not present in this file
            raise ORCAError(self.inputFile)
        #print self.inputFile, self.HOMO_1, self.HOMO, self.LUMO, self.LUMO_1

        
    def calculateTransferIntegral(self, singleChromos):
        epsilon1 = singleChromos[self.chromo1ID].HOMO
        epsilon2 = singleChromos[self.chromo2ID].HOMO
        HOMOSplitting = self.HOMO-self.HOMO_1
        deltaE = epsilon2-epsilon1
        if deltaE**2 > HOMOSplitting**2:
            # print "\n"
            # print self.inputFile
            # print "HOMO Splitting =", HOMOSplitting
            # print "Delta E =", deltaE
            # raw_input("Complex Transfer Integral")
            self.needKoopmans = 1
            return 0
        else:
            self.needKoopmans = 0
        return 0.5*(np.sqrt(np.abs((self.HOMO - self.HOMO_1)**2 - (epsilon2 - epsilon1)**2))) # Positive anyway so np.abs unnecessary
        # return 0.5*(np.sqrt((self.HOMO - self.HOMO_1)**2 - (epsilon2 - epsilon1)**2))


def readCSVData(fileName, mode):
    singleChromoDict = {}
    pairChromoList = []
    document = csv.reader(fileName, delimiter = ',')
    for row in document:
        chromo = chromophore(None, None, singleChromos = None, chromoData = row)
        if mode == 'singles':
            singleChromoDict[int(float(chromo.chromo1ID))] = chromo
        elif mode == 'pairs':
            pairChromoList.append(chromo)
    if mode == 'singles':
        return singleChromoDict
    elif mode == 'pairs':
        return pairChromoList
    
        
def prepareCSVData(singleChromoDict = None, pairChromoList = None):
    singleChromoCSVData = [] # each row = 1 chromo: [ID, x, y, z, HOMO-1, HOMO, LUMO, LUMO+1, chromoLength]
    pairChromoCSVData = [] # each row = 1 pair: [ID1, ID2, HOMO-1, HOMO, LUMO, LUMO+1, Tij]
    if singleChromoDict != None:
        singleChromoKeys = sorted(singleChromoDict.keys())
        for singleChromoKey in singleChromoKeys:
            chromophore = singleChromoDict[singleChromoKey]
            singleChromoCSVData.append([chromophore.chromo1ID, chromophore.position[0], chromophore.position[1], chromophore.position[2], chromophore.HOMO_1, chromophore.HOMO, chromophore.LUMO, chromophore.LUMO_1, chromophore.chromoLength])
    if pairChromoList != None:
        for chromophore in pairChromoList:
            pairChromoCSVData.append([chromophore.chromo1ID, chromophore.chromo2ID, chromophore.HOMO_1, chromophore.HOMO, chromophore.LUMO, chromophore.LUMO_1, chromophore.Tij])
    return np.array(singleChromoCSVData), np.array(pairChromoCSVData)
        
    
def scaleEnergies(singleChromoDict):
    # Shorter chromophores have significantly deeper HOMOs because they are treated as small molecules instead of chain segments.
    # To rectify this, find the average HOMO level for each chromophore length and then map that average to the literature P3HT HOMO
    litHOMO = -5.0 #eV
    # Calculate the average HOMO for each length
    HOMOLengthDict = {}
    for singleChromoKey in singleChromoDict.keys():
        chromophore = singleChromoDict[singleChromoKey]
        if chromophore.chromoLength not in HOMOLengthDict:
            HOMOLengthDict[chromophore.chromoLength] = [chromophore.HOMO]
        else:
            HOMOLengthDict[chromophore.chromoLength].append(chromophore.HOMO)
    # Calculate the change required to map the energy levels to experiment
    deltaHOMO = {}
    for length, HOMO in HOMOLengthDict.iteritems():
        avHOMO = np.average(HOMO)
        deltaHOMO[length] = litHOMO - avHOMO
    # Perform the map
    for singleChromoKey in singleChromoDict.keys():
        chromoLength = singleChromoDict[singleChromoKey].chromoLength
        singleChromoDict[singleChromoKey].HOMO_1 += deltaHOMO[chromoLength]
        singleChromoDict[singleChromoKey].HOMO += deltaHOMO[chromoLength]
        singleChromoDict[singleChromoKey].LUMO += deltaHOMO[chromoLength]
        singleChromoDict[singleChromoKey].LUMO_1 += deltaHOMO[chromoLength]
    # Now shrink down the DoS to 100 meV
    targetstd = 0.1
    HOMOLevels = []
    for singleChromoKey in singleChromoDict.keys():
        HOMOLevels.append(singleChromoDict[singleChromoKey].HOMO)
    mean = np.mean(np.array(HOMOLevels))
    std = np.std(np.array(HOMOLevels))
    for singleChromoKey in singleChromoDict.keys():
        sigma = (singleChromoDict[singleChromoKey].HOMO - mean)/float(std)
        newDeviation = targetstd*sigma
        singleChromoDict[singleChromoKey].HOMO = mean + newDeviation
    newHOMOLevels = []
    for singleChromoKey in singleChromoDict.keys():
        newHOMOLevels.append(singleChromoDict[singleChromoKey].HOMO)
    return singleChromoDict

def recalculateTij(pairChromos, singleChromoDict):
    koopmans = 0
    for chromoPair in pairChromos:
        chromo1ID = chromoPair.chromo1ID
        chromo2ID = chromoPair.chromo2ID
        HOMOSplitting = chromoPair.HOMO-chromoPair.HOMO_1
        deltaE = singleChromoDict[chromo2ID].HOMO - singleChromoDict[chromo1ID].HOMO
        if deltaE**2 > HOMOSplitting**2:
            chromoPair.Tij = 0
            koopmans += 1
        else:
            chromoPair.Tij = 0.5*np.sqrt((HOMOSplitting**2) - (deltaE**2))
    return pairChromos, koopmans
    
    
        
def getChromoID(fileName):
    slashList = helperFunctions.findIndex(fileName, '/')
    if slashList != None:
        fileName = fileName[slashList[-1]+1:]
    return map(int, ('_'+fileName[:-4]).split('_chromo')[1:])


def turnOffSOSCF(inputFile):
    with open(inputFile, 'r') as fileName:
        originalLines = fileName.readlines()
    originalLines[3] = '!ZINDO/S NoSOSCF\n'
    with open(inputFile, 'w+') as fileName:
        fileName.writelines(originalLines)


def reduceTolerance(inputFile):
    with open(inputFile, 'r') as fileName:
        originalLines = fileName.readlines()
    originalLines[3] = '!ZINDO/S NoSOSCF SloppySCF\n'
    with open(inputFile, 'w+') as fileName:
        fileName.writelines(originalLines)


def increaseIterations(inputFile):
    with open(inputFile, 'r') as fileName:
        originalLines = fileName.readlines()
    originalLines.append('\n%scf MaxIter 500 end')
    with open(inputFile, 'w+') as fileName:
        fileName.writelines(originalLines)


def increaseGrid(inputFile):
    with open(inputFile, 'r') as fileName:
        originalLines = fileName.readlines()
    originalLines[3] = '!ZINDO/S SlowConv Grid7 NoFinalGrid\n'
    with open(inputFile, 'w+') as fileName:
        fileName.writelines(originalLines)


def increaseGridNoSOSCF(inputFile):
    with open(inputFile, 'r') as fileName:
        originalLines = fileName.readlines()
    originalLines[3] = '!ZINDO/S SlowConv Grid7 NoFinalGrid NoSOSCF SloppySCF\n'
    originalLines.append('\n%scf MaxIter 500 end')
    with open(inputFile, 'w+') as fileName:
        fileName.writelines(originalLines)

        
def revertORCAFiles(inputFile):
    with open(inputFile, 'r') as fileName:
        originalLines = fileName.readlines()
    originalLines[3] = '! ZINDO/S\n'
    for lineNo in range(len(originalLines)):
        # REMOVE THE SCF ITER
        if "%scf MaxIter" in originalLines[lineNo]:
            originalLines.pop(lineNo)
            break
    with open(inputFile, 'w+') as fileName:
        fileName.writelines(originalLines)


def execute(morphologyFile):
    morphologyName = morphologyFile[helperFunctions.findIndex(morphologyFile,'/')[-1]+1:]
    orcaOutputDir = os.getcwd()+'/outputFiles/'+morphologyName+'/chromophores/outputORCA'
    CSVDir = os.getcwd()+'/outputFiles/'+morphologyName+'/chromophores'
    orcaDir = os.getenv('ORCA_BIN', str(os.getcwd())+'/ORCA')
    orcaPath = orcaDir+'/orca'
    singleCSVPresent = False
    pairCSVPresent = False
    pairCheckDict = {}
    # Load a pickle of the transfer integral dictionaries?
    try:
        with open(CSVDir+'/singles.csv', 'r') as CSVFile:
            singlesCSVPresent = True
            singleChromoDict = readCSVData(CSVFile, 'singles')
    except IOError:
        singleCSVPresent = False
        singleChromoDict = {} # Is a dictionary because we need to be able to look up energy levels quickly
    try:
        with open(CSVDir+'/pairs.csv', 'r') as CSVFile:
            pairCSVPresent = True
            pairChromos = readCSVData(CSVFile, 'pairs')
            for pair in pairChromos:
                pairCheckDict[str([pair.chromo1ID, pair.chromo2ID])] = True
    except IOError:
        pairCSVPresent = False
        pairChromos = []
        
    singleDir = orcaOutputDir+'/single'
    pairDir = orcaOutputDir+'/pair'
    try:
        procIDs = list(np.arange(int(os.environ.get('SLURM_NPROCS'))))
    except (AttributeError, TypeError):
        # Was not loaded using SLURM, so use all physical processors
        procIDs = list(np.arange(mp.cpu_count()))

    failedSinglesDict = {}
    # Load all of the single ORCA chromophores first to get their energy levels
    # (We need these first to calculate Tij for the pairs)
    failedSingleFiles = []
    failedSingleNos = []
    singleOutputs = []
    for fileName in os.listdir(singleDir):
        if '.out' in fileName:
            singleOutputs.append(fileName)
    for fileNo, fileName in enumerate(singleOutputs):
        chromoID = getChromoID(fileName)
        if chromoID[0] in singleChromoDict:
            print "Single chromophore", chromoID, "already present in CSV file. Skipping...\r",
            continue
        print "Examining single chromophore", fileNo+1, "of", str(len(singleOutputs))+"...\r",
        chromo = chromophore(singleDir+'/'+fileName, chromoID)
        if chromo.error == 0:
            singleChromoDict[chromoID[0]] = chromo
        else:
            failedSingleFiles.append(singleDir.replace('outputORCA', 'inputORCA')+'/'+fileName.replace('.out', '.inp'))
            failedSingleNos.append(chromoID[0])
    print "\n"
    print "There were", len(failedSingleNos), "failed single chromophore runs"
    # Modify the energy levels to get the experimental DoS distribution
    singleChromoDict = scaleEnergies(singleChromoDict)
    # Fixed the failed files
    while len(failedSingleFiles) > 0:
        fixedFilesIndices = []
        try:
            print "Calculations completed for single-segments, however there were", len(failedSingleFiles), "errors in calculating the single-segment HOMO levels."
            for failIndex, fileName in enumerate(failedSingleFiles):
                print "Re-examining file:", fileName
                if fileName not in failedSinglesDict:
                    failedSinglesDict[fileName] = 1
                elif failedSinglesDict[fileName] == 18:
                    continue
                else:
                    failedSinglesDict[fileName] += 1
                failedRerunCounter = failedSinglesDict[fileName]
                if failedRerunCounter == 3:
                    # Three lots of reruns without any successes, try to turn off SOSCF
                    print str(fileName)+": Three lots of reruns without any success - turning off SOSCF to see if that helps..."
                    turnOffSOSCF(fileName)
                if failedRerunCounter == 6:
                    # Still no joy - increase the number of SCF iterations and see if convergence was just slow
                    print str(fileName)+": Six lots of reruns without any success - increasing the number of SCF iterations to 500..."
                    increaseIterations(fileName)
                if failedRerunCounter == 9:
                    # Finally, turn down the SCF tolerance
                    print str(fileName)+": Nine lots of reruns without any success - decreasing SCF tolerance (sloppySCF)..."
                    reduceTolerance(fileName)
                if failedRerunCounter == 12:
                    print str(fileName)+": Failed to rerun ORCA 12 times, one final thing that can be done is to change the numerical accuracy..."
                    revertORCAFiles(fileName)
                    increaseGrid(fileName)
                if failedRerunCounter == 15:
                    print str(fileName)+": Failed to rerun ORCA 15 times. Will try high numerical accuracy with no SOSCF as a last-ditch effort..."
                    increaseGridNoSOSCF(fileName)
                if failedRerunCounter == 18:
                    # SERIOUS PROBLEM
                    print str(fileName)+": Failed to rerun ORCA 18 times, even with all the input file tweaks. Examine the geometry - it is most likely unreasonable."
                    print "Reverting "+str(fileName)+" back to its original state..."
                    revertORCAFiles(fileName)
                    fixedFilesIndices.append(failIndex)
                    failedSinglesDict.pop(fileName)
            if len(failedSingleFiles) == 0:
                break
            jobsList = [failedSingleFiles[i:i+(int(np.ceil(len(failedSingleFiles)/len(procIDs))))+1] for i in xrange(0, len(failedSingleFiles), int(np.ceil(len(failedSingleFiles)/float(len(procIDs)))))]
            with open(CSVDir+'/ORCAJobs.pickle', 'w+') as pickleFile:
                pickle.dump(jobsList, pickleFile)
            # Now rerun ORCA
            if len(jobsList) <= len(procIDs):
                procIDs = procIDs[:len(jobsList)]
            runningJobs = []
            for CPURank in procIDs:
                print 'python '+os.getcwd()+'/code/singleCoreRunORCA.py '+os.getcwd()+'/outputFiles/'+morphologyName+' '+str(CPURank)+' &'
                #os.system('python '+os.getcwd()+'/code/singleCoreRunORCA.py '+os.getcwd()+'/outputFiles/'+morphologyName+' '+str(CPURank)+' &')
                runningJobs.append(sp.Popen(['python', str(os.getcwd())+'/code/singleCoreRunORCA.py', str(os.getcwd())+'/outputFiles/'+morphologyName, str(CPURank), '1'])) # The final argument here tells ORCA to ignore the presence of the output file and recalculate
                    # orcaJob = sp.Popen([str(orcaPath), str(fileName)], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
                    # orcaShellOutput = orcaJob.communicate()
                    # helperFunctions.writeToFile(fileName.replace('inputORCA', 'outputORCA').replace('.inp', '.out'), orcaShellOutput[0].split('\n'), mode = 'outputFile')
            # Wait for running jobs to finish
            t0 = T.time()
            exitCodes = [p.wait() for p in runningJobs]
            t1 = T.time()
            # Now check if the file is fixed
            for failIndex, fileName in enumerate(failedSingleFiles):
                outputFileName = fileName.replace('inputORCA', 'outputORCA').replace('.inp', '.out')
                shortFileName = outputFileName[helperFunctions.findIndex(outputFileName, '/')[-1]+1:]
                chromoID = getChromoID(shortFileName)
                chromo = chromophore(outputFileName, chromoID)
                if chromo.error == 0:
                    fixedFilesIndices.append(failIndex)
                    failedSinglesDict.pop(fileName)
                    singleChromoDict[chromoID[0]] = chromo
            # At end of loop, pop all the fixed files
            if len(fixedFilesIndices) > 0:
                for index in sorted(fixedFilesIndices, reverse=True):
                    failedSingleFiles.pop(index)
        except KeyboardInterrupt:
            print "Kill command recieved. Reverting ORCA files..."
            for inputName, failCount in failedSinglesDict.items():
                revertORCAFiles(inputName)
            print "File reversion complete. Terminating..."
            sys.exit(0)
    # Now write the singles.csv file now we have all the data
    singleChromoDict = scaleEnergies(singleChromoDict)
    singleChromoCSVData, emptyVar = prepareCSVData(singleChromoDict, None)
    helperFunctions.writeCSV(CSVDir+'/singles.csv', singleChromoCSVData)

    print "Calculations completed for single-segments with no errors."
    os.system('rm '+str(CSVDir)+'/ORCAJobs.pickle')


    # Now do the pairs
    failedPairsDict = {}
    failedPairFiles = []
    needKoopmans = 0
    pairOutputs = []
    for fileName in os.listdir(pairDir):
        if '.out' in fileName:
            pairOutputs.append(fileName)
    for fileNo, fileName in enumerate(pairOutputs):
        chromoID = getChromoID(fileName)
        if str(chromoID) in pairCheckDict:
            print "Chromophore pair", chromoID, "already present in CSV file. Skipping...\r",
            continue
        print "Examining chromophore pair", fileNo+1, "of", str(len(pairOutputs))+"...\r",
        if (chromoID[0] in failedSingleNos) or (chromoID[1] in failedSingleNos):
            print "\n"
            print "One of", chromoID, "failed in the singles. Skipping..."
            continue
        chromoPair = chromophore(pairDir+'/'+fileName, chromoID, singleChromoDict)
        if chromoPair.error == 0:
            pairChromos.append(chromoPair)
            dataToWrite = [chromoPair.chromo1ID, chromoPair.chromo2ID, chromoPair.HOMO_1, chromoPair.HOMO, chromoPair.LUMO, chromoPair.LUMO_1, chromoPair.Tij]
            helperFunctions.appendCSV(CSVDir+'/pairs.csv', dataToWrite)
            needKoopmans += chromoPair.needKoopmans
        else:
            failedPairFiles.append(pairDir.replace('outputORCA', 'inputORCA')+'/'+fileName.replace('.out', '.inp'))
    # Recalculate the transfer integrals if the csv file is already present (in case scaling the single energies has affected things)
    if pairCSVPresent == True:
        pairChromos, needKoopmans = recalculateTij(pairChromos, singleChromoDict)
        emptyVar, pairCSVData = prepareCSVData(None, pairChromos)
        helperFunctions.writeCSV(CSVDir+'/pairs.csv', pairCSVData)

    print "\n"
    print failedPairFiles
    print needKoopmans, "pairs out of", len(pairChromos), "had DeltaE > HOMO splitting and so need Koopmans approximation to avoid complex Tij"
    print "There were", len(failedPairFiles), "failed pair chromophore runs"
    # Fix the failed files
    while len(failedPairFiles) > 0:
        fixedFilesIndices = []
        try:
            print "Calculations completed for segment pairs, however there were", len(failedPairFiles), "errors in calculating the pair HOMO splitting."
            for failIndex, fileName in enumerate(failedPairFiles):
                if fileName not in failedPairsDict:
                    failedPairsDict[fileName] = 1
                else:
                    failedPairsDict[fileName] += 1
                failedRerunCounter = failedPairsDict[fileName]
                if failedRerunCounter == 3:
                    # Three lots of reruns without any successes, try to turn off SOSCF
                    print str(fileName)+": Three lots of reruns without any success - turning off SOSCF to see if that helps..."
                    turnOffSOSCF(fileName)
                if failedRerunCounter == 6:
                    # Still no joy - increase the number of SCF iterations and see if convergence was just slow
                    print str(fileName)+": Six lots of reruns without any success - increasing the number of SCF iterations to 500..."
                    increaseIterations(fileName)
                if failedRerunCounter == 9:
                    # Finally, turn down the SCF tolerance
                    print str(fileName)+": Nine lots of reruns without any success - decreasing SCF tolerance (sloppySCF)..."
                    reduceTolerance(fileName)
                if failedRerunCounter == 12:
                    print str(fileName)+": Failed to rerun ORCA 12 times, one final thing that can be done is to change the numerical accuracy..."
                    revertORCAFiles(fileName)
                    increaseGrid(fileName)
                if failedRerunCounter == 15:
                    print str(fileName)+": Failed to rerun ORCA 15 times. Will try high numerical accuracy with no SOSCF as a last-ditch effort..."
                    increaseGridNoSOSCF(fileName)
                if failedRerunCounter == 18:
                    # SERIOUS PROBLEM
                    print str(fileName)+": Failed to rerun ORCA 18 times, even with all the input file tweaks. Examine the geometry - it is most likely unreasonable."
                    print "Reverting "+str(fileName)+" back to its original state..."
                    revertORCAFiles(fileName)
                    fixedFilesIndices.append(failIndex)
                    failedPairsDict.pop(fileName)
            # If failed to run 18 times, just continue and skip
            if len(failedPairFiles) == 0:
                break
            jobsList = [failedPairFiles[i:i+(int(np.ceil(len(failedPairFiles)/len(procIDs))))+1] for i in xrange(0, len(failedPairFiles), int(np.ceil(len(failedPairFiles)/float(len(procIDs)))))]
            with open(CSVDir+'/ORCAJobs.pickle', 'w+') as pickleFile:
                pickle.dump(jobsList, pickleFile)
            # NOW RERUN ORCA
            if len(jobsList) <= len(procIDs):
                procIDs = procIDs[:len(jobsList)]
            runningJobs = []
            for CPURank in procIDs:
                print 'python '+os.getcwd()+'/code/singleCoreRunORCA.py '+os.getcwd()+'/outputFiles/'+morphologyName+' '+str(CPURank)+' &'
                runningJobs.append(sp.Popen(['python', str(os.getcwd())+'/code/singleCoreRunORCA.py', str(os.getcwd())+'/outputFiles/'+morphologyName, str(CPURank), '1']))# The final argument here tells ORCA to ignore the presence of the output file and recalculate
            # for failIndex, fileName in enumerate(failedPairFiles):
            #     orcaJob = sp.Popen([str(orcaPath), str(fileName)], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
            #     orcaShellOutput = orcaJob.communicate()
            #     helperFunctions.writeToFile(fileName.replace('inputORCA', 'outputORCA').replace('.inp', '.out'), orcaShellOutput[0].split('\n'), mode = 'outputFile')
            # Sleep for 20 minutes
            exitCodes = [p.wait() for p in runningJobs]
            # Now check if the file is fixed
            for failIndex, fileName in enumerate(failedPairFiles):
                outputFileName = fileName.replace('inputORCA', 'outputORCA').replace('.inp', '.out')
                chromoID = getChromoID(outputFileName)
                chromoPair = chromophore(outputFileName, chromoID, singleChromoDict)
                if chromoPair.error == 0:
                    fixedFilesIndices.append(failIndex)
                    failedPairsDict.pop(fileName)
                    pairChromos.append(chromoPair)
                    # Put in a check here to see why certain attributes don't exist
                    try:
                        dataToWrite = [chromoPair.chromo1ID, chromoPair.chromo2ID, chromoPair.HOMO_1, chromoPair.HOMO, chromoPair.LUMO, chromoPair.LUMO_1, chromoPair.Tij]
                    except AttributeError:
                        print "outputFileName =", outputFileName
                        print "chromoID =", chromoID
                        print "chromoPair attributes =", dir(chromoPair)
                        continue
                    helperFunctions.appendCSV(CSVDir+'/pairs.csv', dataToWrite)
            # At end of loop, pop all the fixed files
            if len(fixedFilesIndices) > 0:
                for index in sorted(fixedFilesIndices, reverse=True):
                    failedPairFiles.pop(index)
        except KeyboardInterrupt:
            print "Kill command recieved. Reverting ORCA files..."
            for inputName, failCount in failedPairsDict.items():
                revertORCAFiles(inputName)
            print "File reversion complete. Terminating..."
            sys.exit(0)

    print "Calculations completed for segment pairs with no errors."
    os.system('rm '+str(CSVDir)+'/ORCAJobs.pickle')


    
    # Now write the CSV outputs for the KMC code
    # singleChromoCSVData, pairChromoCSVData = prepareCSVData(singleChromoDict, pairChromos)
    # helperFunctions.writeCSV(CSVDir+'/singles.csv', singleChromoCSVData)
    # helperFunctions.writeCSV(CSVDir+'/pairs.csv', pairChromoCSVData)


if __name__ == "__main__":
    morphologyFile = sys.argv[1]
    execute(morphologyFile)
