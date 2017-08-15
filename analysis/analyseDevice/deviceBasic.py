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
            # Some problem in the Results section (maybe it didn't finish?)
            # Run the incomplete data parser instead
            fileData = parseDataIncomplete(data)
            if fileData is not None:
                for key, val in fileData.items():
                    dataDict[key].append(val)
    return dataDict


def parseDataIncomplete(data):
    dataDict = {'J': 0, 'V': 0, 'Photo': 0, 'Cathode': 0, 'Anode': 0, 'Diss': 0, 'Rec': 0, 'Ext': 0, 'Iter': 0, 'SimT': 0, 'WallT': 0}
    numberOfExtractions = []
    timeOfExtractions = []
    dataDict['V'] = float(eval(data[1][:-1])[0])
    simT = 0
    wallT = 0
    iterations = 0
    for line in data:
        if 'EVENT' in line:
            if ('Dark Current' and 'Anode' in line):
                dataDict['Anode'] += 1
            elif ('Dark Current' and 'Cathode' in line):
                dataDict['Cathode'] += 1
            elif ('New number of extractions' in line):
                dataDict['Ext'] += 1
                splitLine = line.split(' ')
                extractionNumber = int(splitLine[12])
                extractionTime = float(splitLine[-1][:-2])
                if extractionNumber is not 0:
                    numberOfExtractions.append(extractionNumber)
                    timeOfExtractions.append(extractionTime)
            elif ('Recombination Succeeded' in line):
                dataDict['Rec'] += 1
            elif ('Photoinjetion' in line):
                dataDict['Photo'] += 1
            elif ('Exciton Dissociating' in line):
                dataDict['Diss'] += 1
        elif ('Current runtime' in line):
            wallT = int(line.split(' ')[3][:-2])
            simT = float(line.split(' ')[24][:-3])
            iterations = int(line.split(' ')[20])

    dataDict['SimT'] = simT
    dataDict['WallT'] = wallT
    dataDict['Iter'] = iterations
    if len(numberOfExtractions) != 0:
        gradient, intercept, rVal, pVal, stdErr = scipy.stats.linregress(timeOfExtractions, numberOfExtractions)
        jVal = float((elementaryCharge * gradient) / (deviceArea * 10))
        dataDict['J'] = jVal
        plt.figure()
        plt.scatter(timeOfExtractions, numberOfExtractions)
        plt.xlabel('Time (s)')
        plt.ylabel('Number of Extractions (Arb. U)')
        plt.xlim([0, 2E-4])
        plt.show()
    else:
        return None
    return dataDict


def calculateJ(data, deviceArea):
    numberOfExtractions = [0]
    timeOfExtractions = [0.0]
    for line in data:
        if 'number of extractions' in line:
            splitLine = line.split(' ')
            extractionNumber = int(splitLine[12])
            extractionTime = float(splitLine[-1][:-2])
            if extractionNumber != numberOfExtractions[-1]:
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
    xVals, yVals = zip(*sorted(zip(xVals, yVals)))
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
