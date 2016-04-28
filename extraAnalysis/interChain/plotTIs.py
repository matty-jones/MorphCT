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



def plotHist(yvals, mode, xvals=None):
    plt.figure()
    if mode == 'intraChain':
        plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Intra-Chain TI (eV)')
        fileName = 'intraTij.png'
    elif mode == 'interChain':
        plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Inter-Chain TI (eV)')
        fileName = 'interTij.png'
    if mode == 'intraChainTrim':
        plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Intra-Chain TI (eV)')
        fileName = 'intraTijTrim.png'
    elif mode == 'interChainTrim':
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
    
