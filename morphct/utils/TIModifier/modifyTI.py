import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import csv


def loadCSVs():
    CSVDir = os.getcwd()
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


def writeNewCSV(pairsData):
    CSVDir = os.getcwd()
    with open(CSVDir+'/KoopPairs.csv', 'w+') as pairsFile:
        pairsWriter = csv.writer(pairsFile, delimiter=',')
        for row in pairsData:
            pairsWriter.writerow([str(row[0])+'.0', str(row[1])+'.0', float(row[2]), float(row[3]), float(row[4]), float(row[5]), float(row[6])])
    print "New CSV written to", CSVDir+'/KoopPairs.csv'


def applyKoopmans(singlesData, pairsData):
    newPairData = []
    for i, pair in enumerate(pairsData):
        # Use the HOMO_1 and HOMO levels to calculate the new Tij
        # Tij = 0.5*(HOMO - HOMO-1)
        pairsData[i][6] = 0.5*(pairsData[i][3] - pairsData[i][2])
    writeNewCSV(pairsData)

    
if __name__ == "__main__":
    singleData, pairData = loadCSVs()
    # singleData = {ChromoID: [x, y, z, HOMO-1, HOMO, LUMO, LUMO+1, Length]}
    # pairData = [chromo1, chromo2, HOMO-1, HOMO, LUMO, LUMO+1, TransferIntegral]
    applyKoopmans(singleData, pairData)
