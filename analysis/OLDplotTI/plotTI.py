import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import csv

from scipy.optimize import curve_fit


def loadCSVs(CSVDir):
    singlesData = {}
    pairsData = []
    molIDs = {}
    with open(CSVDir+'/singles.csv', 'r') as singlesFile:
        singlesReader = csv.reader(singlesFile, delimiter=',')
        for row in singlesReader:
            singlesData[int(float(row[0]))] = map(float, row[1:])
    with open(CSVDir+'/pairs.csv', 'r') as pairsFile:
        pairsReader = csv.reader(pairsFile, delimiter=',')
        for row in pairsReader:
            pairsData.append([int(float(row[0])), int(float(row[1]))] + [x for x in map(float, row[2:])])
    with open(CSVDir+'/molIDs.csv', 'r') as molIDsFile:
        molIDsReader = csv.reader(molIDsFile, delimiter=',')
        for row in molIDsReader:
            for chromoID in row[1:]:
                molIDs[int(chromoID)] = int(row[0])
    return singlesData, pairsData, molIDs


def calculateLambdaij(chromoLength):
    # The equation for the internal reorganisation energy was obtained from the data given in
    # Johansson, E and Larsson, S; 2004, Synthetic Metals 144: 183-191.
    # External reorganisation energy obtained from 
    # Liu, T and Cheung, D. L. and Troisi, A; 2011, Phys. Chem. Chem. Phys. 13: 21461-21470
    lambdaExternal = 0.11 # eV
    if chromoLength < 12:
        lambdaInternal = 0.20826 - (chromoLength*0.01196)
    else:
        lambdaInternal = 0.06474
    lambdaeV = lambdaExternal+lambdaInternal
    return lambdaeV


def gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


def gaussFit(data):
    n = len(data)
    mean = np.mean(data)
    std = np.std(data)
    print "\n"
    print "Delta Eij stats: mean =", mean, "std =", std
    hist, binEdges = np.histogram(data, bins=100)
    fitArgs, fitConv = curve_fit(gaussian, binEdges[:-1], hist, p0=[1, mean, std])
    return binEdges, fitArgs


def plotHist(saveDir, yvals, mode, xvals=None, gaussBins=None, fitArgs=None):
    if mode == 'HOMO':
        plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('HOMO Level (eV)')
        fileName = 'HOMODoS.pdf'
    elif mode == 'Bandgap':
        plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Bandgap (eV)')
        fileName = 'Bandgap.pdf'
    elif mode == 'BandgapLength':
        plt.scatter(xvals, yvals)
        plt.ylabel('Bandgap (eV)')
        plt.xlabel('Chromo Length (monomers)')
        fileName = 'BandgapLength.pdf'
    elif mode == 'Splitting':
        plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('HOMO Splitting (ev)')
        fileName = 'HOMOSplit.pdf'
    elif mode == 'TI':
        plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Transfer Integral (eV)')
        fileName = 'TI.pdf'
    elif mode == 'Trimmed':
        plt.hist(yvals, bins = np.arange(0, 1.05, 0.05))
        plt.ylabel('Frequency')
        plt.xlabel(r'$T_{ij}$ (eV)')
        plt.xlim([0, 1.0])
        fileName = 'TITrimmed.pdf'
    elif mode == 'Length':
        plt.scatter(xvals, yvals)
        plt.xlabel('HOMO Level (eV)')
        plt.ylabel('Chromo Length (monomers)')
        fileName = 'HOMOLength.pdf'
    elif mode == 'lambda':
        plt.scatter(xvals, yvals)
        plt.xlabel('Chromo Length (monomers)')
        plt.ylabel('Reorganisation energy (eV)')
        fileName = 'LambdaIJ.pdf'
    elif mode == 'deltaEij':
        n, bins, patches = plt.hist(yvals, 20)
        gaussY = gaussian(gaussBins[:-1], *fitArgs)
        scaleFactor = max(n)/max(gaussY)
        plt.plot(gaussBins[:-1], gaussY*scaleFactor, 'ro:')
        plt.ylabel('Frequency')
        plt.xlabel('Delta Eij (eV)')
        plt.xlim([-1.5, 1.5])
        fileName = 'deltaEij.pdf'
    elif mode == 'averageHOMO':
        plt.scatter(xvals, yvals)
        plt.xlabel('Chromo Length (monomers)')
        plt.ylabel('Average HOMO level (eV)')
        fileName = 'averageHOMO.pdf'
    elif mode == 'intraChain':
        if len(yvals) > 0:
            plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Intra-Chain TI (eV)')
        fileName = 'intraTij.pdf'
    elif mode == 'interChain':
        if len(yvals) > 0:
            plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Inter-Chain TI (eV)')
        fileName = 'interTij.pdf'
    elif mode == 'intraChainTrim':
        if len(yvals) > 0:
            plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Intra-Chain TI (eV)')
        fileName = 'intraTijTrim.pdf'
    elif mode == 'interChainTrim':
        if len(yvals) > 0:
            plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Inter-Chain TI (eV)')
        fileName = 'interTijTrim.pdf'
    elif mode == 'wobbey':
        plt.scatter(xvals, yvals)
        plt.xlabel('Chromo Length (monomers)')
        plt.ylabel('A (eV)')
        fileName = 'averageHOMO.pdf'
    plt.savefig(saveDir+'/'+fileName)
    plt.clf()
    print "Figure saved as", saveDir+"/"+fileName


def findIndex(string, character):
    '''This function returns the locations of an inputted character in an inputted string'''
    index = 0
    locations = []
    while index < len(string):
        if string[index] == character:
            locations.append(index)
        index += 1
    if len(locations) == 0:
        return None
    return locations


if __name__ == "__main__":
    tempDirs = []
    for fileName in os.listdir(os.getcwd()):
        if fileName[0] == 'T':
            tempDirs.append(fileName)
    plt.figure()
    for tempDir in tempDirs:
        CSVDir = os.getcwd()+'/'+tempDir+'/01mon/TIZero'
        singleData, pairData, molIDs = loadCSVs(CSVDir)
        # singleData = {ChromoID: [x, y, z, HOMO-1, HOMO, LUMO, LUMO+1, Length]}
        # pairData = [chromo1, chromo2, HOMO-1, HOMO, LUMO, LUMO+1, TransferIntegral]
        # molIDs = {ChromoID: molID}
        chromoLength = []
        HOMOLevels = []
        bandgap = []
        for chromophore in singleData.values():
            chromoLength.append(chromophore[7])
            HOMOLevels.append(chromophore[4])
            bandgap.append(chromophore[5] - chromophore[4])
        HOMOSplitting = []
        transferIntegrals = []
        trimmedTIs = []
        deltaEij = []
        for chromoPair in pairData:
            HOMOSplitting.append(chromoPair[3]-chromoPair[2])
            transferIntegrals.append(chromoPair[6])
            if chromoPair[6] != 0.0:
                trimmedTIs.append(chromoPair[6])
            deltaEij.append((singleData[chromoPair[0]][4])-(singleData[chromoPair[1]][4]))

        binEdges, fitArgs = gaussFit(deltaEij)

            
        plotHist(CSVDir, HOMOLevels, 'HOMO')
        plotHist(CSVDir, HOMOSplitting, 'Splitting')
        plotHist(CSVDir, transferIntegrals, 'TI')
        plotHist(CSVDir, trimmedTIs, 'Trimmed')
        plotHist(CSVDir, chromoLength, 'Length', xvals=HOMOLevels)
        plotHist(CSVDir, deltaEij, 'deltaEij', gaussBins=binEdges, fitArgs=fitArgs)
        plotHist(CSVDir, bandgap, 'Bandgap')
        plotHist(CSVDir, bandgap, 'BandgapLength', xvals=chromoLength)
        chromoLengths = np.arange(16)
        lambdaIJ = []
        for monomers in chromoLengths:
            lambdaIJ.append(calculateLambdaij(monomers))
        plotHist(CSVDir, lambdaIJ, 'lambda', xvals=chromoLengths)

        HOMODict = {}
        for i, HOMO in enumerate(HOMOLevels):
            length = chromoLength[i]
            if length not in HOMODict:
                HOMODict[chromoLength[i]] = [HOMO]
            else:
                HOMODict[chromoLength[i]].append(HOMO)
        chromophoreLength = []
        averageHOMO = []
        for length, HOMOs in HOMODict.iteritems():
            chromophoreLength.append(length)
            averageHOMO.append(np.average(HOMOs))
        plotHist(CSVDir, averageHOMO, 'averageHOMO', xvals=chromophoreLength)

        interChain = []
        intraChain = []
        interChainTrim = []
        intraChainTrim = []
        for chromoPair in pairData:
            molID = molIDs[chromoPair[0]]
            if molIDs[chromoPair[1]] == molID:
                intraChain.append(chromoPair[6])
                if chromoPair[6] != 0.0:
                    intraChainTrim.append(chromoPair[6])
            else:
                interChain.append(chromoPair[6])
                if chromoPair[6] != 0.0:
                    interChainTrim.append(chromoPair[6])
        plotHist(CSVDir, intraChain, 'intraChain')
        plotHist(CSVDir, interChain, 'interChain')
        plotHist(CSVDir, intraChainTrim, 'intraChainTrim')
        plotHist(CSVDir, interChainTrim, 'interChainTrim')


