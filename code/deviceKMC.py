import os
import sys
import numpy as np
import matplotlib.pyplot as plt
try:
    import mpl_toolkits.mplot3d as p3
except ImportError:
    print("Could not import 3D plotting engine, calling the plotDeviceComponents function will result in an error!")

import helperFunctions


class morphologyMoiety:
    def __init__(self, moietyTypeNumber, molMorphName, parameterDict):
        self.typeNo = moietyTypeNumber
        chromophoreListLocation = parameterDict['outputMorphDir'] + '/' + molMorphName + '/code/' + molMorphName + '.pickle'
        AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, self.parameterDict, self.chromophoreList = helperFunctions.loadPickle(chromophoreListLocation)
        self.carrierType = self.getCarrierType()
        self.TIDist, self.deltaEDist, self.avHopDistance = self.getDistributions(AAMorphologyDict)

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

    def getDistributions(self, AAMorphologyDict):
        boxSize = [AAMorphologyDict['lx'], AAMorphologyDict['ly'], AAMorphologyDict['lz']]
        TIDist = []
        deltaEDist = []
        hopDistanceDist = []
        zValDist = []
        # Iterate over all chromophores in the molecular morphology
        for chromophore in self.chromophoreList:
            for neighbourIndex, neighbourData in enumerate(chromophore.neighbours):
                neighbourID = neighbourData[0]
                neighbourImage = neighbourData[1]
                # Only consider the `forward' hops (to prevent duplication)
                if neighbourID > chromophore.ID:
                    # Ignore hops with zero transfer integral
                    if (chromophore.neighboursTI[neighbourIndex] is not None) and (chromophore.neighboursTI[neighbourIndex] > 0):
                        # Create a transfer integral distribution that we can select from to approximate how carriers move in this moiety
                        TIDist.append(chromophore.neighboursTI[neighbourIndex])
                        # Create a deltaE distribution that we can select from to approximate how carriers move in this moiety
                        deltaEDist.append(chromophore.neighboursDeltaE[neighbourIndex])
                        currentChromoPosn = chromophore.posn
                        neighbourChromoPosn = [self.chromophoreList[neighbourID].posn[axis] + (neighbourImage[axis] * boxSize[axis]) for axis in range(3)]
                        # Create a hopping distance distribution that we can select from to approximate how carriers move in this moiety
                        hopDistanceDist.append(np.linalg.norm(np.array(neighbourChromoPosn) - np.array(currentChromoPosn)))
                        # 
        plotHopDistance(hopDistanceDist)
        return TIDist, deltaEDist, np.average(hopDistanceDist)


def plotHopDistance(distribution):
    plt.figure()
    plt.hist(distribution)
    plt.show()
    exit()




def loadDeviceMorphology(parameterDict):
    deviceDir = parameterDict['inputDeviceDir'] + '/' + parameterDict['deviceMorphology']
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
        moietyDictionary[moietyID] = morphologyMoiety(moietyID, parameterDict['deviceComponents'][moietyID], parameterDict)


    #plotDeviceComponents(deviceArray)


def plotDeviceComponents(deviceArray):
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    colour = ['r', 'g', 'b', 'y', 'k']
    for xVal in range(9):
        for yVal in range(9):
            for zVal in range(9):
                ax.scatter(xVal, yVal, zVal, zdir='z', c= colour[deviceArray[xVal, yVal, zVal]])
    plt.show()
    exit()


def execute(parameterDict):
    # First job will be to load in the device morphology, when I work out what format I want it to be.
    loadDeviceMorphology(parameterDict)
    return


# DEBUG FUNCTION, PLEASE IGNORE
def coldRun(**kwargs):
# DEBUG FUNCTION, PLEASE IGNORE
    parameterDict = {}
    for key, value in kwargs.items():
        parameterDict[key] = value
    execute(parameterDict)


if __name__ == "__main__":
    print("WOBBEY")
