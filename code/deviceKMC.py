import os
import sys
import numpy as np

import helperFunctions


def loadDeviceMorphology(parameterDict):
    deviceDir = parameterDict['inputDeviceDir'] + '/' + parameterDict['deviceMorphology']
    zvals = os.listdir(deviceDir)
    deviceArray = np.zeroes([len(zvals)]*3)
    print deviceArray
    exit()
    for zval in os.listdir(deviceDir):
        pass


def execute(chromophoreList, parameterDict):
    # First job will be to load in the device morphology, when I work out what format I want it to be.
    return


if __name__ == "__main__":
    print("WOBBEY")
