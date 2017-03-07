import numpy as np
import sys
sys.path.append('../../code')
import helperFunctions

if __name__ == "__main__":
    fileName = sys.argv[1]
    # Create the dictionary of the template file using loadMorphologyXML
    morphologyDict = helperFunctions.loadMorphologyXML(fileName)

    morphologyDict = helperFunctions.addUnwrappedPositions(morphologyDict)
    #simDims = [morphologyDict['lx'], morphologyDict['ly'], morphologyDict['lz']]
    #morphologyDict['unwrapped_position'] = []
    #for atomID, posn in enumerate(morphologyDict['position']):
    #    morphologyDict['unwrapped_position'].append(list([morphologyDict['position'][atomID][x] + (simDims[x] * morphologyDict['image'][atomID][x]) for x in range(3)]))

    fullerenesToKeep = np.arange(11249,16060,13)[1:]
    fullerenesToRemove = sorted(list(set(np.arange(11250,16060)) - set(fullerenesToKeep)), reverse=True)
    for constraint in ['position', 'image', 'unwrapped_position', 'mass', 'diameter', 'charge', 'type', 'body']:
        for atomIndex in fullerenesToRemove:
            morphologyDict[constraint].pop(atomIndex)
    morphologyDict['natoms'] -= len(fullerenesToRemove)
    helperFunctions.writeMorphologyXML(morphologyDict, 'EDITTED'+fileName)
