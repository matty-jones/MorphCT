import sys
import os
import csv
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('../../code')
import helperFunctions


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


def createNewCSV(outputFile, data):
    with open(outputFile, 'w+') as csvFile:
        csvWriter = csv.writer(csvFile)
        csvWriter.writerow(['chromoiID', 'chromojID', 'kij', 'chromoiCG', 'chromojCG'])
        for row in data:
            csvWriter.writerow(row)


#def plotkij(inputDir, outputFile, data):
#    intraChainHoppingRates = []
#    interChainHoppingRates = []
#    molIDData = loadCSVData(inputDir+'/molIDs.csv')
#    chromoToMol = {}
#    for moleculeChromos in molIDData:
#        for chromoID in moleculeChromos[1:]:
#            chromoToMol[int(chromoID)] = int(moleculeChromos[0])
#    for hoppingPair in data[1:]:
#        if chromoToMol[hoppingPair[0]] == chromoToMol[hoppingPair[1]]:
#            intraChainHoppingRates.append(hoppingPair[2])
#        else:
#            interChainHoppingRates.append(hoppingPair[2])
#
#    hoppingRates = map(float,list(np.array(data[1:])[:,2]))
#    chromoSeparations = map(float, list(np.array(data[1:])[:,5]))
#
#    plt.figure()
#    plt.hist([intraChainHoppingRates, interChainHoppingRates], bins=np.logspace(1,18,40), stacked=True, color=['red','blue'], label=['Intra-chain', 'Inter-chain'])
#    plt.legend(loc=2, prop={'size': 18})
#    plt.title(outputFile[:-8], fontsize=20)
#    plt.xlim([1,1E18])
#    plt.ylim([0,5000])
#    plt.gca().set_xscale('log')
#    plt.savefig(outputFile.replace('.csv','.pdf'))
#    plt.close()
#
#
#    plt.figure()
#    plt.scatter(chromoSeparations, hoppingRates)
#    plt.title(outputFile[:-8], fontsize=20)
#    plt.xlabel(r'Chromophore Separation / $\AA$')
#    plt.ylabel(r'Hopping Rates / s$^{-1}$')
#    plt.xlim([0, 20])
#    plt.ylim([1, 1E18])
#    plt.gca().set_yscale('log')
#    plt.show()
#    plt.savefig(outputFile.replace('kij.csv','sep.pdf'))
#    plt.close()





if __name__ == "__main__":
    for pickleFile in os.listdir('./inputPickles'):
        outputFile = './outputCSVs/' + pickleFile.replace('.pickle', '_kij.csv')
        AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = helperFunctions.loadPickle('./inputPickles/'+pickleFile)
        CSVData = []
        print("Examining chromophores...")
        for chromo1 in chromophoreList:
            for index, neighbour in enumerate(chromo1.neighbours):
                # Only consider pairs X-Y where Y > X
                chromo2ID = neighbour[0]
                if chromo2ID < chromo1.ID:
                    continue
                Tij = chromo1.neighboursTI[index]
                deltaE = chromo1.neighboursDeltaE[index]
                if (Tij is not None) and (deltaE is not None):
                    Kij = calcHoppingRate(Tij, deltaE)
                    CGIDs1 = chromo1.CGIDs
                    CGIDs2 = chromophoreList[chromo2ID].CGIDs
                    if Kij > 1E5:
                        CSVData.append([str(chromo1.ID), str(chromo2ID), str(Kij), repr(CGIDs1), repr(CGIDs2)])
        print("Writing CSV file...")
        createNewCSV(outputFile, CSVData)
        print("CSV file written to", str(outputFile)+"!")
