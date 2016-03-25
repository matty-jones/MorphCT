import numpy as np
import sys
import helperFunctions
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
import mpl_toolkits.mplot3d.axes3d as p3


def getFunctionalGroups(molecule):
    thioGroups = []
    moleculeEnds = []
    for index, atomType in enumerate(molecule['type']):
        thioGroup = {'index': [], 'type': [], 'mass': [], 'position': []}
        # To be of the form, [IDs], [Types], [Masses], [Positions]
        if atomType == 'S':
            # The S and the preceeding 4 atoms make up the thiophene ring,
            # along with the following H (2 Hs if at the end of the molecule)
            for atomID in range(index-4, index+2):
                thioGroup['index'].append(atomID)
                thioGroup['type'].append(molecule['type'][atomID])
                thioGroup['mass'].append(molecule['mass'][atomID])
                thioGroup['position'].append(molecule['position'][atomID])
            if molecule['type'][index+2] == 'H':
                thioGroup['index'].append(index+2)
                thioGroup['type'].append(molecule['type'][index+2])
                thioGroup['mass'].append(molecule['mass'][index+2])
                thioGroup['position'].append(molecule['position'][index+2])
                moleculeEnds.append(len(thioGroups))
            thioCOM = helperFunctions.calcCOM(thioGroup['position'], thioGroup['mass'])
            thioGroup['COM'] = thioCOM
            # Plane of the thiophene ring is given by the vector between the
            # C1 and C10 (atoms 0 and 3)
            thioPlaneVector = helperFunctions.findAxis(molecule['position'][index-4], molecule['position'][index-1])
            thioGroup['plane'] = thioPlaneVector
            # Normal to plane of thiophene ring is given by the cross product
            # between the C1-C10 vector and the C1-C2 vector
            thioNormalVector = helperFunctions.normaliseVec(np.cross(thioPlaneVector, helperFunctions.findAxis(molecule['position'][index-4], molecule['position'][index-3])))
            thioGroup['normal'] = thioNormalVector
            thioGroups.append(thioGroup)
    return moleculeEnds, thioGroups


def calculateEndToEndDistance(moleculeBackbone, moleculeEnds):
    moleculeEnds = []
    for thioID, thioGroup in enumerate(moleculeBackbone):
        if len(thioGroup['type']) == 7:
            moleculeEnds.append(thioID)
    endToEndDistance = helperFunctions.calculateSeparation(moleculeBackbone[moleculeEnds[0]]['COM'], moleculeBackbone[moleculeEnds[1]]['COM'])
    return endToEndDistance


def calculatePersistanceLength(molecule, molName, outputDir, moleculeBackbone, moleculeEnds):
    distances = []
    cosThetas = []
    for thioRing1 in moleculeBackbone:
        for thioRing2 in moleculeBackbone:
            separation = helperFunctions.calculateSeparation(thioRing1['COM'], thioRing2['COM'])
            distances.append(separation)
            cosThetas.append(np.dot(thioRing1['plane'], thioRing2['plane']))
    distances, cosThetas = helperFunctions.parallelSort(distances, cosThetas)
    startOfCurrentBin = 0
    binWidth = 3.0
    binnedDistances = [[]]
    binnedCosThetas = [[]]
    for i in range(len(distances)):
        if distances[i] > startOfCurrentBin + binWidth:
            binnedDistances.append([])
            binnedCosThetas.append([])
            startOfCurrentBin = distances[i]
        binnedDistances[-1].append(distances[i])
        binnedCosThetas[-1].append(cosThetas[i])
    distanceToPlot = []
    cosThetaToPlot = []
    ybars = []
    for binID in range(len(binnedDistances)):
        distanceToPlot.append(np.average(binnedDistances[binID]))
        cosThetaToPlot.append(np.average(binnedCosThetas[binID]))
        ybars.append(np.std(binnedCosThetas[binID])/np.sqrt(len(binnedCosThetas[binID])))
    plt.figure()
    plt.plot(distanceToPlot, cosThetaToPlot, 'ro')
    plt.errorbar(distanceToPlot, cosThetaToPlot, yerr=ybars, linestyle='None')
    plt.savefig(outputDir+'/molecules/'+molName.replace('.POSCAR', '_PL.png'))
    plt.close()
    return distanceToPlot, cosThetaToPlot


def plotMolecule3D(molecule, moleculeBackbone):
    allAtoms = np.array(molecule['position'])
    COMs = []
    normals = []
    for thioRing in moleculeBackbone:
        normals.append([])
        COMs.append(np.array(thioRing['COM']))
        for i in range(10):
            normals[-1].append(np.array(thioRing['COM']+(thioRing['normal']*i)))
    COMs = np.array(COMs)
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    ax.scatter(allAtoms[:,0], allAtoms[:,1], allAtoms[:,2], s = 20, c = 'g')
    ax.scatter(COMs[:,0], COMs[:,1], COMs[:,2], s = 200, c = 'r')
    for thioRing in normals:
        for coords in thioRing:
            ax.scatter(coords[0], coords[1], coords[2], s = 20, c = 'b')
    ax.set_xlim((-8,8))
    ax.set_ylim((-40,40))
    ax.set_zlim((-15,15))
    plt.show()


def fixAngles(inputAngle):
    if (inputAngle <= np.pi/2.):
            pass
    elif (inputAngle > np.pi/2.) and (inputAngle <= np.pi):
        inputAngle = np.pi - inputAngle
    elif (inputAngle > np.pi) and (inputAngle <= 3*np.pi/2.):
        inputAngle = inputAngle - np.pi
    elif (inputAngle > 3*np.pi/2.) and (inputAngle <= 2*np.pi):
        inputAngle = 2*np.pi - inputAngle
    else:
        print inputAngle
        raise SystemError('ANGLE CALCULATION FAIL')
    return inputAngle

    
def calculateChromophores(molecule, molName, outputDir, moleculeBackbone):
    bendingAngleTolerance = np.pi/6. # The bending angle, beyond which the pi conjugation is said to be broken
    torsionAngleTolerance = np.pi/6. # The torsional angle, beyond which, the pi conjugation is said to be broken
    previousThioPlaneAxis = moleculeBackbone[0]['plane']
    previousThioNormalAxis = moleculeBackbone[0]['normal']
    bendingAngles = []
    torsionAngles = []
    chromophores = [[0]]
    for thioNumberMinus1, thioRing in enumerate(moleculeBackbone[1:]):
        thioNumber = thioNumberMinus1 + 1
        currentThioPlaneAxis = thioRing['plane']
        currentThioNormalAxis = thioRing['normal']
        bendingAngle = fixAngles(np.arccos(np.dot(previousThioPlaneAxis, currentThioPlaneAxis)))
        torsionAngle = fixAngles(np.arccos(np.dot(previousThioNormalAxis, currentThioNormalAxis)))
        if (bendingAngle > bendingAngleTolerance) or (torsionAngle > torsionAngleTolerance):
            # Conjugation is broken, so make a new segment
            chromophores.append([])
        chromophores[-1].append(thioNumber)
        bendingAngles.append(bendingAngle)
        torsionAngles.append(torsionAngle)
        previousThioPlaneAxis = np.copy(currentThioPlaneAxis)
        previousThioNormalAxis = np.copy(currentThioNormalAxis)
    plotHist(bendingAngles, outputDir+'/molecules/'+molName.replace('.POSCAR', '_Bend.png'))
    plotHist(torsionAngles, outputDir+'/molecules/'+molName.replace('.POSCAR', '_Tor.png'))
    return bendingAngles, torsionAngles, chromophores


def plotHist(data, outputFile, bins=20):
    plt.figure()
    plt.hist(data, bins, normed=1)
    plt.savefig(outputFile)
    plt.close()


if __name__ == "__main__":
    morphologyFile = sys.argv[1]
    morphologyName = morphologyFile[helperFunctions.findIndex(morphologyFile,'/')[-1]+1:]
    outputDir = './outputFiles'
    morphologyList = os.listdir(outputDir)
    pickleFound = False
    for allMorphologies in morphologyList:
        if morphologyName in allMorphologies:
            outputDir += '/'+morphologyName
            break
    # for fileName in os.listdir(outputDir+'/morphology'):
    #     if fileName == morphologyName+'.pickle':
    #         pickleLoc = outputDir+'/morphology/'+fileName
    #         pickleFound = True
    # if pickleFound == False:
    #     print "Pickle file not found. Please run morphCT.py again to create the required HOOMD inputs."
    #     exit()
    # print "Pickle found at", str(pickleLoc)+"."
    # print "Loading data..."
    # with open(pickleLoc, 'r') as pickleFile:
    #     (AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, boxSize) = pickle.load(pickleFile)
    # calculateBendingAngles(AAMorphologyDict)







        
    poscarList = os.listdir(outputDir+'/molecules')
    moleculeList = []
    molNames = []
    for poscarFile in poscarList:
        if '.POSCAR' in poscarFile:
            moleculeDict = helperFunctions.loadPoscar(outputDir+'/molecules/'+poscarFile)
            moleculeDict = helperFunctions.addMasses(moleculeDict)
            moleculeList.append(moleculeDict)
            molNames.append(poscarFile)
    morphologyBackboneData = [] # Dictionaries of molecule data for each molecule in the system
    for molNo, molecule in enumerate(moleculeList):
        print "Examining molecule number", molNo, "of", str(len(moleculeList)-1)+"..."
        moleculeEnds, moleculeBackbone = getFunctionalGroups(molecule)
        endToEndDistance = calculateEndToEndDistance(moleculeBackbone, moleculeEnds)
        distances, cosThetas = calculatePersistanceLength(molecule, molNames[molNo], outputDir, moleculeBackbone, moleculeEnds)
        bendingAngles, torsionAngles, chromophores = calculateChromophores(molecule, molNames[molNo], outputDir, moleculeBackbone)
        moleculeDictionary = {'ends': moleculeEnds, 'endToEndDistance': endToEndDistance, 'persistenceLength': [distances, cosThetas], 'bendingAngles': bendingAngles, 'torsionAngles': torsionAngles, 'chromophores': chromophores}
        morphologyBackboneData.append(moleculeDictionary)
    print "Plotting morphology histograms..."
    plotHist([molDict['endToEndDistance'] for molDict in morphologyBackboneData], outputDir+'/morphology/EndToEndDistances.png')
    plotHist([molDict['bendingAngles'] for molDict in morphologyBackboneData], outputDir+'/morphology/BendingAngles.png')
    plotHist([molDict['torsionAngles'] for molDict in morphologyBackboneData], outputDir+'/morphology/TorsionAngles.png')
    print "Average chromophore length =", np.average([len(molDict['chromophores']) for molDict in morphologyBackboneData]) 
