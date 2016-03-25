import numpy as np
import sys
import helperFunctions
import pylab as P
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
    P.figure()
    P.plot(distanceToPlot, cosThetaToPlot, 'ro')
    P.errorbar(distanceToPlot, cosThetaToPlot, yerr=ybars, linestyle='None')
    P.savefig(outputDir+'/molecules/'+molName.replace('.POSCAR', '_PL.png'))
    P.close()
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
    fig = P.figure()
    ax = p3.Axes3D(fig)
    ax.scatter(allAtoms[:,0], allAtoms[:,1], allAtoms[:,2], s = 20, c = 'g')
    ax.scatter(COMs[:,0], COMs[:,1], COMs[:,2], s = 200, c = 'r')
    for thioRing in normals:
        for coords in thioRing:
            ax.scatter(coords[0], coords[1], coords[2], s = 20, c = 'b')
    ax.set_xlim((-8,8))
    ax.set_ylim((-40,40))
    ax.set_zlim((-15,15))
    P.show()

    
def calculateTorsionalAngles(molecule, molName, outputDir, moleculeBackbone):
    referenceAxis = moleculeBackbone[0]['normal']
    torsions = []
    for thioRing in moleculeBackbone[1:]:
        normalAxis = thioRing['normal']
        theta = np.arccos(np.dot(referenceAxis, normalAxis))
        if (theta <= np.pi/2.):
            pass
        elif (theta > np.pi/2.) and (theta <= np.pi):
            theta = np.pi - theta
        elif (theta > np.pi) and (theta <= 3*np.pi/2.):
            theta = theta-np.pi
        elif (theta > 3*np.pi/2.) and (theta <= 2*np.pi):
            theta = 2*np.pi - theta
        else:
            print referenceAxis, normalAxis, theta
            raise SystemError('ANGLE CALCULATION FAIL')
        torsions.append(theta)
        referenceAxis = np.copy(normalAxis)
    P.figure()
    P.hist(torsions, 20, normed=1)
    P.savefig(outputDir+'/molecules/'+molName.replace('.POSCAR', '_Tor.png'))
    P.close()
    return torsions


def plotHist(data, outputFile, bins=20):
    P.figure()
    P.hist(data, bins, normed=1)
    P.savefig(outputFile)
    P.close()


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
    poscarList = os.listdir(outputDir+'/molecules')
    moleculeList = []
    molNames = []
    for poscarFile in poscarList:
        if '.POSCAR' in poscarFile:
            moleculeDict = helperFunctions.loadPoscar(outputDir+'/molecules/'+poscarFile)
            moleculeDict = helperFunctions.addMasses(moleculeDict)
            moleculeList.append(moleculeDict)
            molNames.append(poscarFile)
    endToEndDistances = []
    morphologyTorsions = []
    for molNo, molecule in enumerate(moleculeList):
        print "Examining molecule number", molNo, "of", str(len(moleculeList)-1)+"..."
        moleculeEnds, moleculeBackbone = getFunctionalGroups(molecule)
        endToEndDistances.append(calculateEndToEndDistance(moleculeBackbone, moleculeEnds))
        distances, cosThetas = calculatePersistanceLength(molecule, molNames[molNo], outputDir, moleculeBackbone, moleculeEnds)
        morphologyTorsions += calculateTorsionalAngles(molecule, molNames[molNo], outputDir, moleculeBackbone)
    print "Plotting morphology histograms..."
    plotHist(endToEndDistances, outputDir+'/morphology/EndToEndDistribution.png')
    plotHist(morphologyTorsions, outputDir+'/morphology/TorsionDistribution.png')
