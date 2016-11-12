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


def loadORCAOutput(fileName):
    with open(fileName, 'r') as orcaFile:
        dataFile = orcaFile.readlines()
    recordMOData = False
    orbitalData = []
    for line in dataFile:
        if 'ORBITAL ENERGIES' in line:
            # Next line begins the MO data
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
    for i in range(len(orbitalData)):
        if orbitalData[i][1] == 0:
            # This line is the first unoccupied orbital - i.e. LUMO
            LUMO = orbitalData[i][3]
            HOMO = orbitalData[i-1][3]
            HOMO_1 = orbitalData[i-2][3]
            LUMO_1 = orbitalData[i+1][3]
            # Don't need any other orbitals
            break
    if recordMOData == False:
        # Molecular orbital data not present in this file
        raise ORCAError(self.inputFile)
    return HOMO_1, HOMO, LUMO, LUMO_1


def rerunFails(failedChromoFiles):
    # Need a clever function that works out which ones have failed and updates the ones that have been fixed
    pass



def updateChromophoreList(chromophoreList, parameterDict):
    orcaOutputDir = parameterDict['outputDir'] + '/' + parameterDict['morphology'][:-4] + '/chromophores/outputORCA'
    # NOTE: This can possibly be done by recursively iterating through the neighbourlist of each chromophore, but I
    # imagine Python will whinge about the levels of recursion, so for now I'll just go through every chromophore twice.
    # Firstly, set the energy levels for each single chromophore, rerunning them if they fail.
    failedSingleChromos = {}
    for chromophore in chromophoreList:
        fileName = '/single/%04d.out' % (chromophore.ID)
        # Update the chromophores in the chromophoreList with their energyLevels
        try:
            chromophore.HOMO_1, chromophore.HOMO, chromophore.LUMO, chromophore.LUMO_1 = loadORCAOutput(fileName)
        except ORCAError:
            failedSingleChromos[fileName] = 1
    # Rerun any failed ORCA jobs
    while len(failedSingleChromos) > 0:
        failedSingleChromos = rerunFails(failedSingleChromos)

    # Now that all the single chromophore energy levels are done, iterate through again and check the neighbours,
    # rerunning the pair file if it failed (which it won't have done because all my chromophores are delicious now).


    # Delete the orca output files when they're done based on the parameterDict variable (MAKE THIS)




def execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList)
    updateChromophoreList(chromophoreList, parameterDict)


    
    # Now write the CSV outputs for the KMC code
    # singleChromoCSVData, pairChromoCSVData = prepareCSVData(singleChromoDict, pairChromos)
    # helperFunctions.writeCSV(CSVDir+'/singles.csv', singleChromoCSVData)
    # helperFunctions.writeCSV(CSVDir+'/pairs.csv', pairChromoCSVData)

if __name__ == "__main__":
    try:
        pickleFile = sys.argv[1]
    except:
        print "Please specify the pickle file to load to continue the pipeline from this point."
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList = helperFunctions.loadPickle(pickleFile)
    execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList)










































































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
