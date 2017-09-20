import sys
import numpy as np
sys.path.append('../../../code')
import helperFunctions

if __name__ == "__main__":
    fileName = './emptyMorph.xml'
    morphDict = helperFunctions.loadMorphologyXML(fileName)
    morphDict['lx'] = 100.0
    morphDict['ly'] = 100.0
    morphDict['lz'] = 100.0
    for xVal in np.arange(-45.0, 50.0, 10.0):
        for yVal in np.arange(-45.0, 50.0, 10.0):
            for zVal in np.arange(-45.0, 50.0, 10.0):
                morphDict['position'].append([xVal, yVal, zVal])
                morphDict['mass'].append(1.0)
                morphDict['image'].append([0, 0, 0])
                morphDict['body'].append(len(morphDict['position']) - 1)
                if zVal > 0:
                    morphDict['type'].append('D')
                else:
                    morphDict['type'].append('A')
                morphDict['diameter'].append(1.0)
                morphDict['charge'].append(0.0)
    morphDict['natoms'] = len(morphDict['position'])
    helperFunctions.writeMorphologyXML(morphDict, './mixedCrystalBilayer.xml')
