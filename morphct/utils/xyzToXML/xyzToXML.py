import os
import sys
import numpy as np
sys.path.append('../../code')
import helperFunctions


if __name__ == "__main__":
    filesToConvert = {}
    fileNames = os.listdir('./')
    for fileName in fileNames:
        if fileName[-4:] == '.xml':
            filesToConvert[fileName] = []
    for fileName in fileNames:
        if fileName[-4:] == '.xyz':
            filesToConvert[fileName[:helperFunctions.findIndex(fileName, '_')[-1]] + '.xml'].append(fileName)
    for key, vals in filesToConvert.iteritems():
        baseDict = helperFunctions.loadMorphologyXML('./' + key)
        for xyz in vals:
            xyzData = map(list, np.loadtxt('./' + xyz, skiprows = 2, dtype = {'names': ['type', 'x', 'y', 'z'], 'formats': ['S1', 'f4', 'f4', 'f4']}))
            natoms = len(xyzData)
            for atomNo, atom in enumerate(xyzData):
                image = [0, 0, 0]
                posn = atom[1:]
                for axis, key in enumerate(['lx', 'ly', 'lz']):
                    while posn[axis] > baseDict[key]/2.0:
                        posn[axis] -= baseDict[key]
                        image[axis] += 1
                    while posn[axis] < -baseDict[key]/2.0:
                        posn[axis] += baseDict[key]
                        image[axis] -= 1
                baseDict['position'][atomNo] = posn
                baseDict['image'][atomNo] = image
            baseDict['natoms'] = natoms
            helperFunctions.writeMorphologyXML(baseDict, xyz.replace('.xyz', '.xml'))
