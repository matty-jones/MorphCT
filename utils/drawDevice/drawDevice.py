import sys
import os
import numpy as np
sys.path.append('../../code')
import helperFunctions
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as p3


class morphologyMoiety:
    def __init__(self, molMorphName):
        chromophoreListLocation = './' + molMorphName + '/code/' + molMorphName + '.pickle'
        self.AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, self.parameterDict, self.chromophoreList = helperFunctions.loadPickle(chromophoreListLocation)
        self.carrierType = self.getCarrierType()
        # Now add the occupation data to the chromophoreLists so that we can prevent double occupation in the simulations
        # The occupied property is a list that contains the device moiety coordinates where the chromophore is occupied.
        for index, chromophore in enumerate(self.chromophoreList):
            chromophore.occupied = []

    def getCarrierType(self):
        speciesPresent = []
        for chromophore in self.chromophoreList:
            speciesPresent.append(chromophore.species)
        if len(set(speciesPresent)) == 1:
            if speciesPresent[0] == 'Donor':
                return 'hole'
            elif speciesPresent[0] == 'Acceptor':
                return 'electron'
            else:
                print("Error in chromophore:")
                for key, val in chromophore.__dict__.items():
                    print(key, "=", val)
                raise SystemError("Chromophore species is neither Donor nor Acceptor")
        else:
            return 'both'


class chromophoreDataContainer:
    # A helper class that contains all of the chromophore data for ease of access from anywhere
    def __init__(self, deviceArray, moietyDictionary, wrapxy):
        self.deviceArray = deviceArray
        self.moietyDictionary = moietyDictionary

    def returnChromophoreList(self, deviceMoietyType):
        return self.moietyDictionary[deviceMoietyType].chromophoreList



def loadDeviceMorphology(deviceName, deviceComponents):
    deviceDir = './' + deviceName
    ySlices = os.listdir(deviceDir)
    # Initialize the array of the correct size (assumes cubic morphology)
    deviceArray = np.zeros([len(ySlices)]*3, dtype=int)
    for yVal, fileName in enumerate(ySlices):
        # Load the ySlice as-presented in the input files
        ySlice = np.loadtxt(deviceDir + '/' + fileName, dtype=int)
        # The z-origin is at the top, and we need it at the bottom, so turn the array upside down
        ySlice = np.flipud(ySlice)
        # Now populate the array
        for zVal, zRow in enumerate(ySlice):
            for xVal, datum in enumerate(zRow):
                deviceArray[xVal, yVal, zVal] = datum
    moietyDictionary = {}
    for moietyID in np.unique(deviceArray):
        moietyDictionary[moietyID] = morphologyMoiety(deviceComponents[moietyID])
    return deviceArray, moietyDictionary


def plotDevice(deviceArray, chromophoreDict):
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    deviceShape = deviceArray.shape
    colour = {'Donor': 'r', 'Acceptor': 'b'}
    maxZVal = 0
    minZVal = 0
    for deviceX in range(deviceShape[0]):
        for deviceY in range(deviceShape[1]):
            for deviceZ in range(deviceShape[2]):
                for chromophore in chromophoreDict[deviceArray[deviceX, deviceY, deviceZ]]:
                    [xVal, yVal, zVal] = np.array(chromophore.posn) + 100.0 * (np.array([deviceX, deviceY, deviceZ]) - np.array([1, 1, 1]))
                    if zVal > maxZVal:
                        maxZVal = zVal
                    elif zVal < minZVal:
                        minZVal = zVal
                    species = chromophore.species
                    ax.scatter(xVal, yVal, zVal, zdir='z', c=colour[species])
    print(maxZVal)
    print(minZVal)
    for deviceX in range(deviceShape[0]):
        for deviceY in range(deviceShape[1]):
            [xVal, yVal, zVal] = [0, 0, maxZVal] + 100.0 * (np.array([deviceX, deviceY, 0]) - np.array([1, 1, 0]))
            # ANODE IS GREEN
            ax.scatter(xVal, yVal, zVal, zdir='z', c='g')
    for deviceX in range(deviceShape[0]):
        for deviceY in range(deviceShape[1]):
            [xVal, yVal, zVal] = [0, 0, minZVal] + 100.0 * (np.array([deviceX, deviceY, 0]) - np.array([1, 1, 0]))
            # CATHODE IS BLACK
            ax.scatter(xVal, yVal, zVal, zdir='z', c='k')
    plt.show()


if __name__ == "__main__":
    deviceName = 'testCrystalBilayer'
    deviceComponents = {0: 'donorCrystal', 1: 'acceptorCrystal', 2: 'mixedCrystalBilayer'}
    deviceArray, moietyDictionary = loadDeviceMorphology(deviceName, deviceComponents)
    chromophoreData = chromophoreDataContainer(deviceArray, moietyDictionary, False)
    chromophoreDict = {}
    for deviceMoietyType in deviceComponents.keys():
        chromophoreDict[deviceMoietyType] = chromophoreData.returnChromophoreList(deviceMoietyType)
    plotDevice(deviceArray, chromophoreDict)
