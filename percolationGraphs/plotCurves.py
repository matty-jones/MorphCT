import os
import csv
import matplotlib.pyplot as plt


class morphology:
    def __init__(self, csvFile, row):
        self.csvFile = csvFile
        self.temperature = float(row[0])
        self.largestCluster = int(row[1])
        self.clusterQuantity = int(row[2])
        self.clusterList = map(int, row[3][1:-1].split(', '))


def getClusterDistData(effectiveTemperatureData):
    effTempXVals = []
    clusterSizeYVals = []
    clusterQuantYVals = []
    for effTemp, morphology in effectiveTemperatureData.iteritems():
        

        
def importMobility():
    mobData = {}
    with open('./mobility.csv', 'r') as csvFile:
        csvReader = csv.reader(csvFile)
        for row in csvReader:
            temperature = row[0]
            mobility = row[1]
            mobData[temperature] = mobility
    return mobData

        
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


def plotXY(xvals, yvals, mode, morphologyName = None):
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(xvals, yvals, c='#0000a5')
    if mode == "effTempSize":
        ax.set_xlabel(r"Effective Temperature (Arb. U.)")
        ax.set_ylabel(r"Biggest Cluster Size")
        plotName = "./"+morphologyName+"_effTempSize.pdf"
    elif mode == "effTempQuant":
        ax.set_xlabel(r"Effective Temperature (Arb. U.)")
        ax.set_ylabel(r"Cluster Quantity")
        plotName = "./"+morphologyName+"_effTempQuant.pdf"
    elif mode == "sizeMob":
        ax.set_xlabel(r"Biggest Cluster Size")
        ax.set_ylabel(r"Mobility (cm$^{2}$ V$^{-1}$ s$^{-1}$)")
        plotName = "./"+morphologyName+"_sizeMob.pdf"
    elif mode == "quantMob":
        ax.set_xlabel(r"Cluster Quantity")
        ax.set_ylabel(r"Mobility (cm$^{2}$ V$^{-1}$ s$^{-1}$)")
        plotName = "./"+morphologyName+"_quantMob.pdf"
    fig.savefig(plotName)
    print "Figure saved as", plotName
    plt.clf()




if __name__ == "__main__":
    mobData = importMobility()
    effectiveTempDirs = []
    for dirName in os.listdir(os.getcwd()):
        if 'K' in dirName:
            effectiveTempDirs.append(os.getcwd()+'/'+dirName)
    effectiveTemperatureData = {}
    for effectiveTempDir in effectiveTempDirs:
        effectiveTemperature = int(effectiveTempDir[findIndex(effectiveTempDir, '/')[-1]+1:-1])
        CSVFileName = effectiveTempDir+'/clusterData.csv'
        with open(CSVFileName, 'r') as csvFile:
            csvReader = csv.reader(csvFile)
            for row in csvReader:
                effectiveTemperatureData[effectiveTemperature] = morphology(CSVFileName, row)
        
    # Plot out the variation in Cluster Size and Cluster Quantity for each morphology as a function of effective temp
    clusterSizeXVals, clusterQuantXVals, effectiveTempXVals = getClusterDistData(effectiveTemperatureData)

    # For each effective temperature, plot Cluster Size and Quantity against mobility value to see if there's any sort of trend at all
