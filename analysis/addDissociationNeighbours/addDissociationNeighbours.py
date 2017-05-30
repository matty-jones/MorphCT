import numpy as np
import os
import sys
sys.path.append('../../code')

import helperFunctions

if __name__ == "__main__":
    pickleFiles = []
    for fileName in os.listdir('./'):
        if 'pickle' in fileName:
            pickleFiles.append(fileName)
    for pickleFile in pickleFiles:
        print("Treating", pickleFile)
        AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = helperFunctions.loadPickle(pickleFile)
        simDims = [[-AAMorphologyDict['lx']/2.0, AAMorphologyDict['lx']/2.0], [-AAMorphologyDict['ly']/2.0, AAMorphologyDict['ly']/2.0], [-AAMorphologyDict['lx']/2.0, AAMorphologyDict['lz']/2.0]]
        for chromo1 in chromophoreList:
            for chromo2 in chromophoreList:
                if chromo1.ID == chromo2.ID:
                    continue
                deltaPosn = chromo2.posn - chromo1.posn
                relativeImageOfChromo2 = [0, 0, 0]
                for axis in range(3):
                    halfBoxLength = (simDims[axis][1] - simDims[axis][0]) / 2.0
                    while deltaPosn[axis] > halfBoxLength:
                        deltaPosn[axis] -= simDims[axis][1] - simDims[axis][0]
                        relativeImageOfChromo2[axis] -= 1
                    while deltaPosn[axis] < - halfBoxLength:
                        deltaPosn[axis] += simDims[axis][1] - simDims[axis][0]
                        relativeImageOfChromo2[axis] += 1
                separation = np.linalg.norm(deltaPosn)
                if separation <= parameterDict['maximumHopDistance']:
                    try:
                        chromo1NeighbourIDs = [neighbourData[0] for neighbourData in chromo1.dissociationNeighbours]
                    except AttributeError:
                        chromo1.dissociationNeighbours = []
                    try:
                        chromo2NeighbourIDs = [neighbourData[0] for neighbourData in chromo2.dissociationNeighbours]
                    except AttributeError:
                        chromo2.dissociationNeighbours = []
                    if (chromo1.species != chromo2.species):
                        if chromo2.ID not in chromo1NeighbourIDs:
                            chromo1.dissociationNeighbours.append([chromo2.ID, relativeImageOfChromo2])
                        if chromo1.ID not in chromo2NeighbourIDs:
                            chromo2.dissociationNeighbours.append([chromo1.ID, list(-np.array(relativeImageOfChromo2))])
        helperFunctions.writePickle((AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList), pickleFile)


