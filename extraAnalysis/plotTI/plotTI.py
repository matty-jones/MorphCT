import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import csv


def loadCSVs():
    CSVDir = os.getcwd()
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


def plotHist(yvals, mode, xvals=None):
    plt.figure()
    if mode == 'HOMO':
        plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('HOMO Level (eV)')
        fileName = 'HOMODoS.png'
    elif mode == 'Bandgap':
        plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Bandgap (eV)')
        fileName = 'Bandgap.png'
    elif mode == 'BandgapLength':
        plt.scatter(xvals, yvals)
        plt.ylabel('Bandgap (eV)')
        plt.xlabel('Chromo Length (monomers)')
        fileName = 'BandgapLength.png'
    elif mode == 'Splitting':
        plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('HOMO Splitting (ev)')
        fileName = 'HOMOSplit.png'
    elif mode == 'TI':
        plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Transfer Integral (eV)')
        fileName = 'TI.png'
    elif mode == 'Trimmed':
        plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Non-Zero Transfer Integral (eV)')
        fileName = 'TITrimmed.png'
    elif mode == 'Length':
        plt.scatter(xvals, yvals)
        plt.xlabel('HOMO Level (eV)')
        plt.ylabel('Chromo Length (monomers)')
        fileName = 'HOMOLength.png'
    elif mode == 'lambda':
        plt.scatter(xvals, yvals)
        plt.xlabel('Chromo Length (monomers)')
        plt.ylabel('Reorganisation energy (eV)')
        fileName = 'LambdaIJ.png'
    elif mode == 'deltaEij':
        plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Delta Eij (eV)')
        fileName = 'deltaEij.png'
    elif mode == 'averageHOMO':
        plt.scatter(xvals, yvals)
        plt.xlabel('Chromo Length (monomers)')
        plt.ylabel('Average HOMO level (eV)')
        fileName = 'averageHOMO.png'
    elif mode == 'intraChain':
        if len(yvals) > 0:
            plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Intra-Chain TI (eV)')
        fileName = 'intraTij.png'
    elif mode == 'interChain':
        if len(yvals) > 0:
            plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Inter-Chain TI (eV)')
        fileName = 'interTij.png'
    elif mode == 'intraChainTrim':
        if len(yvals) > 0:
            plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Intra-Chain TI (eV)')
        fileName = 'intraTijTrim.png'
    elif mode == 'interChainTrim':
        if len(yvals) > 0:
            plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Inter-Chain TI (eV)')
        fileName = 'interTijTrim.png'
    elif mode == 'wobbey':
        plt.scatter(xvals, yvals)
        plt.xlabel('Chromo Length (monomers)')
        plt.ylabel('A (eV)')
        fileName = 'averageHOMO.png'
    plt.savefig('./'+fileName)
    print "Figure saved as ./"+fileName


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
    singleData, pairData, molIDs = loadCSVs()
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
    plotHist(HOMOLevels, 'HOMO')
    plotHist(HOMOSplitting, 'Splitting')
    plotHist(transferIntegrals, 'TI')
    plotHist(trimmedTIs, 'Trimmed')
    plotHist(chromoLength, 'Length', HOMOLevels)
    plotHist(deltaEij, 'deltaEij')
    plotHist(bandgap, 'Bandgap')
    plotHist(bandgap, 'BandgapLength', chromoLength)
    chromoLengths = np.arange(16)
    lambdaIJ = []
    for monomers in chromoLengths:
        lambdaIJ.append(calculateLambdaij(monomers))
    plotHist(lambdaIJ, 'lambda', chromoLengths)



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
    plotHist(averageHOMO, 'averageHOMO', chromophoreLength)

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
    plotHist(intraChain, 'intraChain')
    plotHist(interChain, 'interChain')
    plotHist(intraChainTrim, 'intraChainTrim')
    plotHist(interChainTrim, 'interChainTrim')
    

