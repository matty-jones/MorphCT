import sys
import os
import matplotlib.pyplot as plt
import scipy.stats

elementaryCharge = 1.60217662E-19


def loadDataFiles(directory):
    fileList = os.listdir(directory + '/KMC/')
    filesToAnalyse = []
    for fileName in fileList:
        if fileName[-4:] == '.log':
            filesToAnalyse.append(directory + '/KMC/' + fileName)
    return filesToAnalyse


def parseData(dataFileList, deviceArea):
    dataDict = {'J': [], 'V': [], 'Photo': [], 'Cathode': [], 'Anode': [], 'Diss': [], 'Rec': [], 'Ext': [], 'Iter': [], 'SimT': [], 'WallT': []}
    for fileName in dataFileList:
        print(fileName)
        with open(fileName, 'r') as fileHandle:
            data = fileHandle.readlines()
            quickCheck(data)
            JVal = calculateJ(data, deviceArea)  # in (C/s) / m^2
            JVal /= 10.0  # in mA/cm^{2}
        try:
            timeData = data[-7].split(' ')
            dataDict['Iter'].append(int(timeData[3]))
            dataDict['SimT'].append(float(timeData[7][:-1]))
            dataDict['WallT'].append(float(timeData[9]))
            dataDict['V'].append(float(data[-8].split(' ')[-1][:-1]))
            dataDict['J'].append(JVal)
            dataDict['Photo'].append(int(data[-6].split(' ')[-1][:-1]))
            dataDict['Cathode'].append(int(data[-5].split(' ')[-1][:-1]))
            dataDict['Anode'].append(int(data[-4].split(' ')[-1][:-1]))
            dataDict['Diss'].append(int(data[-3].split(' ')[-1][:-1]))
            dataDict['Rec'].append(int(data[-2].split(' ')[-1][:-1]))
            dataDict['Ext'].append(int(data[-1].split(' ')[-1][:-1]))
        except ValueError:
            continue
    return dataDict


def calculateJ(data, deviceArea):
    numberOfExtractions = []
    timeOfExtractions = []
    for line in data:
        if 'number of extractions' in line:
            splitLine = line.split(' ')
            extractionNumber = int(splitLine[12])
            extractionTime = float(splitLine[-1][:-2])
            if extractionNumber is not 0:
                numberOfExtractions.append(extractionNumber)
                timeOfExtractions.append(extractionTime)
    #plt.figure()
    #plt.scatter(timeOfExtractions, numberOfExtractions)
    #plt.xlabel('Time (s)')
    #plt.ylabel('Number of Extractions (Arb. U)')
    #plt.xlim([0, 2E-4])
    #plt.show()
    #exit()
    gradient, intercept, rVal, pVal, stdErr = scipy.stats.linregress(timeOfExtractions, numberOfExtractions)
    return (elementaryCharge * gradient) / deviceArea


def quickCheck(data):
    darkCurrentInjs = {'Anode': [], 'Cathode': []}
    for line in data:
        if 'EVENT: Dark Current' in line:
            darkCurrentData = line.split(' ')
            darkCurrentInjs[darkCurrentData[6]].append(int(darkCurrentData[15][:-1]))


def plotData(xLabel, xVals, yLabel, yVals, fileName, mode='scatter'):
    plt.figure()
    if mode == 'scatter':
        plt.scatter(xVals, yVals)
    else:
        plt.plot(xVals, yVals)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.savefig(fileName)
    print("Plot of", yLabel, "against", xLabel, "saved to", fileName)


if __name__ == "__main__":
    deviceArea = (3 * 7.14E-9)**2
    print("Using a device area of", str(deviceArea) + ". Make sure this is correct for the system that is being studied!")
    deviceDirectory = sys.argv[1]
    dataFiles = loadDataFiles(deviceDirectory)
    dataDict = parseData(dataFiles, deviceArea)
    plotData('Voltage (V)', dataDict['V'], 'J (mA / cm^{2})', dataDict['J'], deviceDirectory + '/Figures/JV.pdf', mode='line')
