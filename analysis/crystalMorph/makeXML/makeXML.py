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
    newAtomList = []
    for xVal in np.arange(-45.0, 50.0, 10.0):
        for yVal in np.arange(-45.0, 50.0, 10.0):
            for zVal in np.arange(-45.0, 50.0, 10.0):
                newAtom = {}
                newAtom['position'] = [xVal, yVal, zVal]
                newAtom['mass'] = 1.0
                newAtom['image'] = [0, 0, 0]
                newAtom['body'] = -1
                if zVal > 0:
                    newAtom['type'] = 'D'
                else:
                    newAtom['type'] = 'A'
                newAtom['diameter'] = 1.0
                newAtom['charge'] = 0.0
                newAtomList.append(newAtom)
                ## In order to calculate the transfer integrals out of the morphology moiety, we also need to include a layer of the `wrong material' that can be picked up by the periodic boundary calculation. Otherwise carriers will permanently be trapped in the mixed bilayer.
                #if zVal == 45.0:
                #    morphDict['position'].append([xVal, yVal, zVal])
                #    morphDict['mass'].append(1.0)
                #    morphDict['image'].append([0, 0, 0])
                #    morphDict['body'].append(len(morphDict['position']) - 1)
                #    morphDict['type'].append('A')
                #    morphDict['diameter'].append(1.0)
                #    morphDict['charge'].append(0.0)
                #elif zVal == -45.0:
                #    morphDict['position'].append([xVal, yVal, zVal])
                #    morphDict['mass'].append(1.0)
                #    morphDict['image'].append([0, 0, 0])
                #    morphDict['body'].append(len(morphDict['position']) - 1)
                #    morphDict['type'].append('D')
                #    morphDict['diameter'].append(1.0)
                #    morphDict['charge'].append(0.0)
    # Now add in the periodic bits
    # Line of Donor just outside the morphology along the top (-ve Z with [0, 0, 1] image)
    for xVal in np.arange(-45.0, 50.0, 10.0):
        for yVal in np.arange(-45.0, 50.0, 10.0):
            zVal = -45.0
            newAtom = {}
            newAtom['position'] = [xVal, yVal, zVal]
            newAtom['mass'] = 1.0
            newAtom['image'] = [0, 0, 1]
            newAtom['body'] = -1
            newAtom['type'] = 'D'
            newAtom['diameter'] = 1.0
            newAtom['charge'] = 0.0
            newAtomList.append(newAtom)
    # Line of Acceptor just outside the morphology along the bottom (+ve Z with [0, 0, -1] image)
    for xVal in np.arange(-45.0, 50.0, 10.0):
        for yVal in np.arange(-45.0, 50.0, 10.0):
            zVal = 45.0
            newAtom = {}
            newAtom['position'] = [xVal, yVal, zVal]
            newAtom['mass'] = 1.0
            newAtom['image'] = [0, 0, -1]
            newAtom['body'] = -1
            newAtom['type'] = 'A'
            newAtom['diameter'] = 1.0
            newAtom['charge'] = 0.0
            newAtomList.append(newAtom)
    # Now populate the morphDict
    # Donor first
    for atom in newAtomList:
        if atom['type'] == 'D':
            for key, val in atom.items():
                morphDict[key].append(val)
    # And now Acceptor
    for atom in newAtomList:
        if atom['type'] == 'A':
            for key, val in atom.items():
                morphDict[key].append(val)
    # Finally, update the rigid bodies
    for index, body in enumerate(morphDict['body']):
        morphDict['body'][index] = index
    morphDict['natoms'] = len(morphDict['position'])
    helperFunctions.writeMorphologyXML(morphDict, './perfectBilayer.xml')
