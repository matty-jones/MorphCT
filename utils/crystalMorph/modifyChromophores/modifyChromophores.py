import sys
import numpy as np
sys.path.append('../../../code')
import helperFunctions

if __name__ == "__main__":
    pickleFile = sys.argv[1]
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = helperFunctions.loadPickle(pickleFile)
    for chromophore in chromophoreList:
        chromophore.HOMO = -5.3
        chromophore.LUMO = -3.9
        chromo1Posn = chromophore.posn
        for index, [neighbourID, neighbourImage] in enumerate(chromophore.neighbours):
            # THE FOLLOWING CODE IS USEFUL IF WE ARE USING THE VORONOI CUTOFF
            chromophore.neighboursTI[index] = 1.0
            # THE FOLLOWING CODE IS USEFUL IF WE ARE USING THE NEIGHBOUR CUTOFF
            #chromo2Posn = np.array(chromophoreList[neighbourID].posn) + (np.array(neighbourImage) * np.array([AAMorphologyDict['lx'], AAMorphologyDict['ly'], AAMorphologyDict['lz']]))
            #separation = helperFunctions.calculateSeparation(chromo1Posn, chromo2Posn)
            ## The following equation results in a TI of 1 eV for the nearest neighbours and
            ## 0.36 eV for the next nearest neighbours (at a separation of root 2)
            #chromophore.neighboursTI[index] = np.e ** (1 - ((separation / 10.0)**2))
            chromophore.neighboursDeltaE[index] = 0.0
    helperFunctions.writePickle((AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList), pickleFile)
