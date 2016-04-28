import numpy as np
import sys
import helperFunctions
import csv


class ORCAError(Exception):
    def __init__(self, fileName):
        self.string = "No molecular orbital data present for "+str(fileName)
    def __str__(self):
        return self.string
        
class chromophore:
    def __init__(self, inputFile, chromoID, singleChromos = None):
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
        deltaE = np.absolute(epsilon1-epsilon2)
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


def prepareCSVData(singleChromoDict, pairChromoList):
    singleChromoCSVData = [] # each row = 1 chromo: [ID, x, y, z, HOMO-1, HOMO, LUMO, LUMO+1]
    pairChromoCSVData = [] # each row = 1 pair: [ID1, ID2, HOMO-1, HOMO, LUMO, LUMO+1, Tij]
    singleChromoKeys = sorted(singleChromoDict.keys())
    for singleChromoKey in singleChromoKeys:
        chromophore = singleChromoDict[singleChromoKey]
        singleChromoCSVData.append([chromophore.chromo1ID, chromophore.position[0], chromophore.position[1], chromophore.position[2], chromophore.HOMO_1, chromophore.HOMO, chromophore.LUMO, chromophore.LUMO_1, chromophore.chromoLength])
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
    return singleChromoDict
    
        
def getChromoID(fileName):
    return map(int, ('_'+fileName[:-4]).split('_chromo')[1:])



def execute(morphologyFile):
    morphologyName = morphologyFile[helperFunctions.findIndex(morphologyFile,'/')[-1]+1:]
    orcaOutputDir = os.getcwd()+'/outputFiles/'+morphologyName+'/chromophores/outputORCA'
    CSVDir = os.getcwd()+'/outputFiles/'+morphologyName+'/chromophores'
    singleCSVPresent = False
    pairCSVPresent = False
    # Load a pickle of the transfer integral dictionaries?
    try:
        with open(CSVDir+'/singles.csv', 'r'):
            singlesCSVPresent = True
    except IOError:
        singleCSVPresent = False
    try:
        with open(CSVDir+'/pairs.csv', 'r'):
            pairCSVPresent = True
    except IOError:
        pairCSVPresent = False

    if (singleCSVPresent == True) and (pairCSVPresent == True):
        print "Data files already present for this morphology. Exitting..."
        return
        
    singleDir = orcaOutputDir+'/single/'
    pairDir = orcaOutputDir+'/pair/'

    # Load all of the single ORCA chromophores first to get their energy levels
    # (We need these first to calculate Tij for the pairs)
    singleChromoDict = {} # Is a dictionary because we need to be able to look up energy levels quickly
    failedSingleFiles = []
    failedSingleNos = []
    singleOutputs = []
    for fileName in os.listdir(singleDir):
        if '.out' in fileName:
            singleOutputs.append(fileName)
    for fileNo, fileName in enumerate(singleOutputs):
        print "Examining single chromophore", fileNo+1, "of", str(len(singleOutputs))+"...\r",
        chromoID = getChromoID(fileName)
        chromo = chromophore(singleDir+fileName, chromoID)
        if chromo.error == 0:
            singleChromoDict[chromoID[0]] = chromo
        else:
            failedSingleFiles.append(fileName)
            failedSingleNos.append(chromoID[0])
    print "\n"
    print failedSingleNos
    print "There were", len(failedSingleNos), "failed single chromophore runs"
    # Modify the energy levels to get the experimental DoS distribution
    for i in singleChromoDict.keys():
        print singleChromoDict[i].HOMO
        break
    singleChromoDict = scaleEnergies(singleChromoDict)
    for i in singleChromoDict.keys():
        print singleChromoDict[i].HOMO
        break
    exit()
    # Now do the pairs
    pairChromos = []
    failedPairFiles = []
    needKoopmans = 0
    pairOutputs = []
    for fileName in os.listdir(pairDir):
        if '.out' in fileName:
            pairOutputs.append(fileName)
    for fileNo, fileName in enumerate(pairOutputs):
        chromoID = getChromoID(fileName)
        print "Examining chromophore pair", fileNo+1, "of", str(len(pairOutputs))+"...\r",
        if (chromoID[0] in failedSingleNos) or (chromoID[1] in failedSingleNos):
            print "\n"
            print "One of", chromoID, "failed in the singles. Skipping..."
            continue
        chromoPair = chromophore(pairDir+fileName, chromoID, singleChromoDict)
        if chromoPair.error == 0:
            pairChromos.append(chromoPair)
            if chromoPair.needKoopmans == 1:
                pass
            needKoopmans += chromoPair.needKoopmans
        else:
            failedPairFiles.append(fileName)
            
    print "\n"
    print failedPairFiles
    print needKoopmans, "pairs out of", len(pairChromos), "had DeltaE > HOMO splitting and so need Koopmans approximation to avoid complex Tij"
    print "There were", len(failedPairFiles), "failed pair chromophore runs"
    # Now write the CSV outputs for the KMC code
    singleChromoCSVData, pairChromoCSVData = prepareCSVData(singleChromoDict, pairChromos)
    helperFunctions.writeCSV(CSVDir+'/singles.csv', singleChromoCSVData)
    helperFunctions.writeCSV(CSVDir+'/pairs.csv', pairChromoCSVData)


if __name__ == "__main__":
    morphologyFile = sys.argv[1]
    execute(morphologyFile)
