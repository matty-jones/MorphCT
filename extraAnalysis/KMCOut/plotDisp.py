import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import csv


def loadCSVs():
    CSVDir = os.getcwd()
    CSVList = []
    completeCSVData = []
    for fileName in os.listdir(CSVDir):
        if ".csv" in fileName:
            CSVList.append(CSVDir+'/'+fileName)
    for CSVFileName in CSVList:
        completeCSVData.append([])
        with open(CSVFileName, 'r') as CSVFile:
            CSVData = csv.reader(CSVFile, delimiter=',')
            for row in CSVData:
                completeCSVData[-1].append(map(float, row))
    return completeCSVData
        
    
def plotHist(CSVFile):
    # CSVArray = np.array(CSVFile)
    # MSD = np.average(CSVArray[:,1]**2)
    # meanTime = np.average(CSVArray[:,3])
    disps = []
    squaredDisps = []
    times = []
    for carrierNo, carrierData in enumerate(CSVFile):
        if carrierData[2] != 1: #Skip single hops
            disps.append(carrierData[1])
            squaredDisps.append(carrierData[1]**2)
            times.append(carrierData[3])
    MSD = np.average(squaredDisps)
    meanTime = np.average(times)
    plt.figure()
    plt.hist(disps, 20)
    plt.ylabel('Frequency')
    plt.xlabel('Displacement (m)')
    fileName = 'disp_%.2E.png' % (meanTime)
    plt.savefig('./'+fileName)
    print "Figure saved as ./"+fileName
    return MSD, meanTime


def plotMSD(times, MSDs):
    plt.figure()
    plt.semilogx(times, MSDs)
    plt.xlabel('Time (s)')
    plt.ylabel('MSD (m^{2})')
    fileName = 'MSD.png'
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
    completeCSVData = loadCSVs()
    times = []
    MSDs = []
    for CSVFile in completeCSVData:
        MSD, meanTime = plotHist(CSVFile)
        times.append(meanTime)
        MSDs.append(MSD)
    plotMSD(times, MSDs)
