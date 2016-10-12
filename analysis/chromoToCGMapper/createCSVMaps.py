import sys
import os
import csv
import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt

def loadCSVData(CSVLoc, mapToInt=False):
    csvLines = []
    with open(CSVLoc, 'r') as csvFile:
        csvReader = csv.reader(csvFile)
        for row in csvReader:
            if mapToInt is True:
                csvLines.append(list(map(int,row)))
            else:
                csvLines.append(row)
    return csvLines


def calcHoppingRate(Tij, deltaEij):
    elementaryCharge = 1.60217657E-19 # C
    lambdaij = 0.3063  # eV
    Tij *= elementaryCharge
    deltaEij *= elementaryCharge
    lambdaij *= elementaryCharge
    kB = 1.3806488E-23 # m^{2} kg s^{-2} K^{-1}
    hbar = 1.05457173E-34 # m^{2} kg s^{-1}
    temperature = 290 # K
    kij = ((2*np.pi)/hbar)*(Tij**2)*np.sqrt(1.0/(4*lambdaij*np.pi*kB*temperature))*np.exp(-((deltaEij+lambdaij)**2)/(4*lambdaij*kB*temperature))
    return kij


def treatData(singlesDict, pairsData, chromoDict):
    headerLine = ['chromoiID', 'chromojID', 'kij', 'chromoiCG', 'chromojCG', 'separation']
    linesToWrite = []
    for hoppingPair in pairsData:
        chromoiID = int(hoppingPair[0])
        chromojID = int(hoppingPair[1])
        deltaEij = singlesDict[chromojID][5] - singlesDict[chromoiID][5]  # [5] == HOMO level
        Tij = float(hoppingPair[6])
        kij = calcHoppingRate(Tij, deltaEij)
        if kij == 0.0:
            continue
        chromoiCG = repr(chromoDict[chromoiID]['CGIDs'])
        chromojCG = repr(chromoDict[chromojID]['CGIDs'])
        chromoiPosn = np.array(singlesDict[int(hoppingPair[0])][1:4])
        chromojPosn = np.array(singlesDict[int(hoppingPair[1])][1:4])
        separation = np.linalg.norm(chromojPosn - chromoiPosn)
        # INCLUDE PERIODIC BOUNDARY CONDITIONS!!!
        boxSize = 23.99707984924
        if separation > boxSize/2.0:
            separation = boxSize - separation
        linesToWrite.append([chromoiID, chromojID, kij, chromoiCG, chromojCG, separation])
    linesToWrite.sort(key = lambda x: int(x[1]))
    linesToWrite.sort(key = lambda x: int(x[0]))
    return [headerLine] + linesToWrite


def createNewCSV(outputFile, data):
    with open(outputFile, 'w+') as csvFile:
        csvWriter = csv.writer(csvFile)
        for row in data:
            csvWriter.writerow(row)


def plotkij(inputDir, outputFile, data):
    intraChainHoppingRates = []
    interChainHoppingRates = []
    molIDData = loadCSVData(inputDir+'/molIDs.csv')
    chromoToMol = {}
    for moleculeChromos in molIDData:
        for chromoID in moleculeChromos[1:]:
            chromoToMol[int(chromoID)] = int(moleculeChromos[0])
    for hoppingPair in data[1:]:
        if chromoToMol[hoppingPair[0]] == chromoToMol[hoppingPair[1]]:
            intraChainHoppingRates.append(hoppingPair[2])
        else:
            interChainHoppingRates.append(hoppingPair[2])

    hoppingRates = map(float,list(np.array(data[1:])[:,2]))
    chromoSeparations = map(float, list(np.array(data[1:])[:,5]))

    plt.figure()
    plt.hist([intraChainHoppingRates, interChainHoppingRates], bins=np.logspace(1,18,40), stacked=True, color=['red','blue'], label=['Intra-chain', 'Inter-chain'])
    plt.legend(loc=2, prop={'size': 18})
    plt.title(outputFile[:-8], fontsize=20)
    plt.xlim([1,1E18])
    plt.ylim([0,5000])
    plt.gca().set_xscale('log')
    plt.savefig(outputFile.replace('.csv','.pdf'))
    plt.close()


    plt.figure()
    plt.scatter(chromoSeparations, hoppingRates)
    plt.title(outputFile[:-8], fontsize=20)
    plt.xlabel(r'Chromophore Separation / $\AA$')
    plt.ylabel(r'Hopping Rates / s$^{-1}$')
    plt.xlim([0, 20])
    plt.ylim([1, 1E18])
    plt.gca().set_yscale('log')
    plt.show()
    plt.savefig(outputFile.replace('kij.csv','sep.pdf'))
    plt.close()





if __name__ == "__main__":
    # inputDir = sys.argv[1]
    # outputFile = sys.argv[2]
    for inputDir in ['T1.5', 'T1.75', 'T2.0', 'T2.25', 'T2.5']:
        outputFile = 'p1-L15-f0.0-P0.1-'+inputDir+'-e0.5_kij.csv'
        with open(inputDir+'/chromophores.pickle') as pickleFile:
            print "Loading pickle..."
            fullChromoDict = pickle.load(pickleFile)
        chromoDict = {}
        for chromoID in fullChromoDict:
            if fullChromoDict[chromoID]['periodic'] is False:
                chromoDict[chromoID] = fullChromoDict[chromoID]
        print "Loading single CSV data..."
        singlesData = loadCSVData(inputDir+'/singles.csv')
        singlesDict = {}
        for row in singlesData:
            singlesDict[int(float(row[0]))] = map(float, row[1:])
        print "Loading pair CSV data..."
        pairsData = loadCSVData(inputDir+'/pairs.csv')
        print "Updating csv data..."
        outputData = treatData(singlesDict, pairsData, chromoDict)
        print "Writing output to", str(outputFile)+"..."
        createNewCSV(outputFile, outputData)
        print "Plotting histogram of hopping rates..."
        plotkij(inputDir, outputFile, outputData)
        print "All done!"
