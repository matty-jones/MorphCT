import numpy as np
import sys
import os
import copy
import matplotlib.pyplot as plt
sys.path.append('../../code/')
import helperFunctions
try:
    import mpl_toolkits.mplot3d.axes3d as p3
except ImportError:
    print "Could not import 3D plotting engine, calling the plotMolecule3D function will result in an error!"
    pass


def loadXMLs(directoryName):
    XMLList = []
    for fileName in os.listdir(directoryName):
        if (fileName[-4:] == '.xml') and ('relaxed_' not in fileName):
            CGMorphology = directoryName+'/'+fileName
            for secondFileName in os.listdir(directoryName):
                if (secondFileName[-4:] == '.xml') and ('relaxed_' in secondFileName) and (fileName[3:-4] in secondFileName):
                    AAMorphology = directoryName+'/'+secondFileName
                    XMLList.append([CGMorphology, AAMorphology])
    morphologyList = []
    namesList = []
    for XMLPair in XMLList:
        morphologyList.append([])
        hyphenLocs = helperFunctions.findIndex(XMLPair[0], '-')
        namesList.append(XMLPair[0][hyphenLocs[-2]+1:hyphenLocs[-1]])
        for XMLFile in XMLPair:
            morphologyList[-1].append(helperFunctions.loadMorphologyXML(XMLFile))
    return morphologyList, namesList

def coarseGrain(AAMorph, CGMorph):
    residue = 0
    while True:
        if (residue*377) + 376 > len(AAMorph['unwrapped_position'])-1:
            break
        # New monomer
        monomer = 0
        previousMonomer = None
        while True:
            if monomer == 15:
                break
            thioRingDict = {'unwrapped_position':[], 'mass':[]}
            alk1Dict = {'unwrapped_position':[], 'mass':[]}
            alk2Dict = {'unwrapped_position':[], 'mass':[]}
            # ThioRing = (residue*377)+(monomer : monomer+4)
            for key in ['unwrapped_position', 'mass']:
                thioRingDict[key] = AAMorph[key][(residue*377)+(monomer*25):(residue*377)+(monomer*25)+5]+[AAMorph[key][(residue*377)+(monomer*25)+24]]
                alk1Dict[key] = AAMorph[key][(residue*377)+(monomer*25)+5:(residue*377)+(monomer*25)+9]+AAMorph[key][(residue*377)+(monomer*25)+18:(residue*377)+(monomer*25)+23]
                alk2Dict[key] = AAMorph[key][(residue*377)+(monomer*25)+8:(residue*377)+(monomer*25)+18]
            thioCOM = helperFunctions.calcCOM(thioRingDict['unwrapped_position'], thioRingDict['mass'])
            CGMorph = appendDictionary(CGMorph, thioCOM, 'A')
            alk1COM = helperFunctions.calcCOM(alk1Dict['unwrapped_position'], alk1Dict['mass'])
            CGMorph = appendDictionary(CGMorph, alk1COM, 'B')
            alk2COM = helperFunctions.calcCOM(alk2Dict['unwrapped_position'], alk2Dict['mass'])
            CGMorph = appendDictionary(CGMorph, alk2COM, 'C')
            # Now add bonds
            if previousMonomer != None:
                CGMorph['bond'].append(['bondA', len(CGMorph['type'])-6, len(CGMorph['type'])-3])
            CGMorph['bond'].append(['bondB', len(CGMorph['type'])-3, len(CGMorph['type'])-2])
            CGMorph['bond'].append(['bondC', len(CGMorph['type'])-2, len(CGMorph['type'])-1])
            previousMonomer = True
            monomer += 1
        residue += 1
    return CGMorph


def appendDictionary(CGMorph, position, atomType):
    CGMorph['position'].append(list(position))
    CGMorph['image'].append([0, 0, 0])
    CGMorph['velocity'].append([0.0, 0.0, 0.0])
    CGMorph['mass'].append(1)
    CGMorph['diameter'].append(1)
    CGMorph['type'].append(atomType)
    CGMorph['body'].append(-1)
    CGMorph['charge'].append(0)
    return CGMorph


def rearrangeDictionary(unorderedDictionary, orderedDictionary, currentIndex, newIndex):
    for key in unorderedDictionary.keys():
        if key in ['lx', 'ly', 'lz', 'natoms', 'bond', 'angle', 'dihedral', 'improper']:
            continue
        # print key, orderedDictionary[key][newIndex], unorderedDictionary[key][currentIndex]
        # print currentIndex, newIndex, len(orderedDictionary[key]), len(unorderedDictionary[key])
        else:
            orderedDictionary[key][newIndex] = copy.deepcopy(unorderedDictionary[key][currentIndex])
    return orderedDictionary


def reIndexBonds(morphology, indexDict):
    for key in ['bond', 'angle', 'dihedral', 'improper']:
        for index, propertyValue in enumerate(morphology[key]):
            for elementNo, element in enumerate(propertyValue):
                if elementNo == 0:
                    continue
                morphology[key][index][elementNo] = indexDict[str(morphology[key][index][elementNo])]
    return morphology



def reorderMorphology(unorderedCGMorph):
    # The initialCGMorph is arranged by molecule such that:
    # Every 45 entries is a new molecule
    # All the thios are the first 15
    # All the alk1s are the second 15
    # All the alk2s are the third 15

    # At this point, unorderedCGMorph is done as 15x Thio-Alk1-Alk2.
    # Therefore, for comparison, it must be reordered.

    orderedCGMorph = copy.deepcopy(unorderedCGMorph)


    # Order list to say when each value comes
    orderList = []
    indexDict = {}
    for moleculeNo in range(250):
        indices = np.arange((moleculeNo*45), ((moleculeNo+1)*45), 3)
        orderList += list(indices) + list(indices+1) + list(indices+2)

    # Now rearrange all the dictionary elements
    for orderedIndex, unorderedIndex in enumerate(orderList):
        orderedCGMorph = rearrangeDictionary(unorderedCGMorph, orderedCGMorph, unorderedIndex, orderedIndex)
        indexDict[str(unorderedIndex)] = orderedIndex

    # Finally, update all the bonds etc.
    print "Refactoring bonds, angles and dihedrals..."
    #orderedCGMorph = reIndexBonds(orderedCGMorph, indexDict)
        
    # Return the new dictionary
    return orderedCGMorph
        
            
def calcDistances(initialMorph, finalMorph):
    distances = []
    for index, pos in enumerate(initialMorph['unwrapped_position']):
        if finalMorph['type'][index] != 'A':
            continue
        separationVector = np.array(finalMorph['position'][index]) - np.array(pos)
        for axis in range(len(separationVector)):
            while separationVector[axis] >= finalMorph['lx']/2.0:
                separationVector[axis] -= finalMorph['lx']
            while separationVector[axis] <= -finalMorph['lx']/2.0:
                separationVector[axis] += finalMorph['lx']
        distance = np.linalg.norm(separationVector)
        if distance > 15:
            print index
            raw_input('PAUSE...')
        distances.append(distance)
    return np.average(distances), distances


def plotHist(distanceDistribution, name):
    plt.hist(distanceDistribution, bins=50)
    plt.xlim([0,5])
    plt.ylim([0, 1000])
    plt.savefig('./distanceDist_'+name+'.pdf')


if __name__ == "__main__":
    plt.figure()
    morphologyList, namesList = loadXMLs(os.getcwd())
    for morphNo, morphologyPair in enumerate(morphologyList):
        print "\nCoarse-graining", namesList[morphNo], "..."
        initialCGMorph = morphologyPair[0]
        finalAAMorph = morphologyPair[1]
        initialCGMorph = helperFunctions.addUnwrappedPositions(initialCGMorph)
        finalAAMorph = helperFunctions.addUnwrappedPositions(finalAAMorph)
        finalCGMorph = {'position':[], 'image':[], 'velocity':[], 'mass':[], 'diameter':[], 'type':[], 'body':[], 'bond':[], 'angle':[], 'dihedral':[], 'improper':[], 'charge':[], 'lx':0, 'ly':0, 'lz':0, 'natoms':0}
        finalCGMorph = coarseGrain(finalAAMorph, finalCGMorph)
        finalCGMorph['lx'] = finalAAMorph['lx']
        finalCGMorph['ly'] = finalAAMorph['ly']
        finalCGMorph['lz'] = finalAAMorph['lz']
        finalCGMorph['natoms'] = len(finalCGMorph['type'])
        CGFileName = './outputCGMorphs/CG_'+namesList[morphNo]+'.xml'
        print "Reordering morphology for consistency..."
        finalCGMorph = reorderMorphology(finalCGMorph)
        print "Scaling morphology to match initial CG..."
        scaleFactor = initialCGMorph['lx']/finalCGMorph['lx']
        finalCGMorph = helperFunctions.scale(finalCGMorph, scaleFactor)
        helperFunctions.writeMorphologyXML(finalCGMorph, CGFileName)
        print "Now determining the separation between before and after atoms..."
        averageDistances, distanceDistribution = calcDistances(initialCGMorph, finalCGMorph)
        print averageDistances
        plt.clf()
        plotHist(distanceDistribution, namesList[morphNo])
    plt.close()
        
