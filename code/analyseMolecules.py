import numpy as np
import sys
import helperFunctions
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
import mpl_toolkits.mplot3d.axes3d as p3


def getFunctionalGroups(molecule, CGtoAAIDs, CGBonds):
    monomerData = []
    for thioID in sorted(CGtoAAIDs.keys()):
        if CGtoAAIDs[thioID][0] == 'thio':
            currentMonomer = [CGtoAAIDs[thioID][1]]
            for bond in CGBonds:
                if (bond[0] == 'bondB') and (bond[1] == thioID):
                    alk1ID = bond[2]
                    currentMonomer.append(CGtoAAIDs[alk1ID][1])
                elif (bond[0] == 'bondB') and (bond[2] == thioID):
                    alk1ID = bond[1]
                    currentMonomer.append(CGtoAAIDs[alk1ID][1])
            for bond in CGBonds:
                if (bond[0] == 'bondC') and (bond[1] == alk1ID):
                    alk2ID = bond[2]
                    currentMonomer.append(CGtoAAIDs[alk2ID][1])
                elif (bond[0] == 'bondC') and (bond[2] == alk1ID):
                    alk2ID = bond[1]
                    currentMonomer.append(CGtoAAIDs[alk2ID][1])
            monomerData.append(currentMonomer)
        else:
            continue
    moleculeEnds = []
    for monomerNo, monomer in enumerate(monomerData):
        if len(monomer[0]) == 7:
            moleculeEnds.append(monomerNo)


    print moleculeEnds
    exit()




    
    print "USE BONDS TO FIND FUNCTIONAL GROUPS"
    exit()
    tempMolAtoms = []
    bondDict = {}
    for bond in molecule['bond']:
        if bond[1] not in bondDict:
            bondDict[bond[1]] = [bond[2]]
        else:
            bondDict[bond[1]].append(bond[2])
        if bond[2] not in bondDict:
            bondDict[bond[2]] = [bond[1]]
        else:
            bondDict[bond[2]].append(bond[1])
    for atomID in range(len(molecule['type'])):
        tempMolAtoms.append([atomID, molecule['type'][atomID], molecule['mass'][atomID], molecule['position'][atomID], sorted(bondDict[atomID])])
    for atom in tempMolAtoms:
        # Start with the first one
        exit()
    
    
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
    hydrogensPerCarbon = {}
    bondTypes = []

    

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
    planes = []
    for thioRing in moleculeBackbone:
        normals.append([])
        planes.append([])
        COMs.append(np.array(thioRing['COM']))
        for i in range(10):
            normals[-1].append(np.array(thioRing['COM']+(np.array(thioRing['normal'])*i/5.)))
            planes[-1].append(np.array(thioRing['COM']+(np.array(thioRing['plane'])*i/5.)))
    COMs = np.array(COMs)
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    #ax.scatter(allAtoms[:,0], allAtoms[:,1], allAtoms[:,2], s = 20, c = 'g')
    ax.scatter(COMs[:,0], COMs[:,1], COMs[:,2], s = 50, c = 'r')
    for thioRing in normals:
        for coords in thioRing:
            ax.scatter(coords[0], coords[1], coords[2], s = 20, c = 'b')
    for thioPlane in planes:
        for coords in thioPlane:
            ax.scatter(coords[0], coords[1], coords[2], s = 20, c = 'g')
    #ax.set_xlim((-8,8))
    ax.set_xlim((-40,40))
    ax.set_ylim((-40,40))
    #ax.set_zlim((-15,15))
    ax.set_zlim((-40,40))
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
    bendingAngleTolerance = np.pi/3. # The bending angle, beyond which the pi conjugation is said to be broken
    torsionAngleTolerance = 2*np.pi/5. # The torsional angle, beyond which, the pi conjugation is said to be broken
    previousThioPlaneAxis = moleculeBackbone[0]['plane']
    previousThioNormalAxis = moleculeBackbone[0]['normal']
    bendingAngles = []
    torsionAngles = []
    chromophores = [[0]]
    conjugationBroken = False
    torsionOverThreshold = 0
    bendingOverThreshold = 0
    for thioNumberMinus1, thioRing in enumerate(moleculeBackbone[1:]):
        thioNumber = thioNumberMinus1 + 1
        currentThioPlaneAxis = thioRing['plane']
        currentThioNormalAxis = thioRing['normal']
        bendingAngle = fixAngles(np.arccos(np.dot(previousThioPlaneAxis, currentThioPlaneAxis)))
        torsionAngle = fixAngles(np.arccos(np.dot(previousThioNormalAxis, currentThioNormalAxis))) # Always flip one of the axes because we have 100% regioregular head-to-tail
        # Create a new segment if conjugation is broken
        if (torsionAngle > torsionAngleTolerance):
            print "Torsion has broken conjugation between", thioNumber-1, "and", thioNumber, ":", torsionAngle*180/np.pi, ">", torsionAngleTolerance*180/np.pi
            torsionOverThreshold += 1
            conjugationBroken = True
        if (bendingAngle > bendingAngleTolerance):
            print "Bending has broken conjugation between", thioNumber-1, "and", thioNumber, ":", bendingAngle*180/np.pi, ">", bendingAngleTolerance*180/np.pi
            bendingOverThreshold += 1
            conjugationBroken = True
        if conjugationBroken == True:
            chromophores.append([])
            conjugationBroken = False
        chromophores[-1].append(thioNumber)
        bendingAngles.append(bendingAngle)
        torsionAngles.append(torsionAngle)
        previousThioPlaneAxis = np.copy(currentThioPlaneAxis)
        previousThioNormalAxis = np.copy(currentThioNormalAxis)
    print "The torsional angle was over the threshold ("+str(torsionAngleTolerance)+")", torsionOverThreshold, "times."
    print "The bending angle was over the threshold ("+str(bendingAngleTolerance)+")", bendingOverThreshold, "times."
    #plotHist(bendingAngles, outputDir+'/molecules/'+molName.replace('.POSCAR', '_Bend.png'), angle=True)
    #plotHist(torsionAngles, outputDir+'/molecules/'+molName.replace('.POSCAR', '_Tor.png'), angle=True)
    return bendingAngles, torsionAngles, chromophores


def plotHist(data, outputFile, bins=20, angle=False):
    if angle == True:
        data=np.array(data)*180/np.pi
        bins = np.arange(0,90,5)
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
    for fileName in os.listdir(outputDir+'/morphology'):
        if fileName == morphologyName+'.pickle':
            pickleLoc = outputDir+'/morphology/'+fileName
            pickleFound = True
    if pickleFound == False:
        print "Pickle file not found. Please run morphCT.py again to create the required HOOMD inputs."
        exit()
    print "Pickle found at", str(pickleLoc)+"."
    print "Loading data..."
    with open(pickleLoc, 'r') as pickleFile:
        (AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, boxSize) = pickle.load(pickleFile)

    CGBonds = []
    for CGBond in CGMoleculeDict['bond']:
        if CGBond not in CGBonds:
            CGBonds.append(CGBond)

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
        moleculeEnds, moleculeBackbone = getFunctionalGroups(molecule, CGtoAAIDs[0], CGBonds)
        # plotMolecule3D(molecule, moleculeBackbone)
        # exit()
        endToEndDistance = calculateEndToEndDistance(moleculeBackbone, moleculeEnds)
        distances, cosThetas = calculatePersistanceLength(molecule, molNames[molNo], outputDir, moleculeBackbone, moleculeEnds)
        bendingAngles, torsionAngles, chromophores = calculateChromophores(molecule, molNames[molNo], outputDir, moleculeBackbone)
        moleculeDictionary = {'ends': moleculeEnds, 'endToEndDistance': endToEndDistance, 'persistenceLength': [distances, cosThetas], 'bendingAngles': bendingAngles, 'torsionAngles': torsionAngles, 'chromophores': chromophores}
        morphologyBackboneData.append(moleculeDictionary)
    print "Plotting morphology histograms..."
    endToEndDistance = []
    bendingAngles = []
    torsionAngles = []
    for molDict in morphologyBackboneData:
        endToEndDistance.append(molDict['endToEndDistance'])
        bendingAngles += molDict['bendingAngles']
        torsionAngles += molDict['torsionAngles']
    plotHist(endToEndDistance, outputDir+'/morphology/EndToEndDistances.png')
    plotHist(bendingAngles, outputDir+'/morphology/BendingAngles.png', angle=True)
    plotHist(torsionAngles, outputDir+'/morphology/TorsionAngles.png', angle=True)
    segmentLengths = []
    for molDict in morphologyBackboneData:
        for segment in molDict['chromophores']:
            segmentLengths.append(len(segment))
    print "Average chromophore length in monomers =", np.average(segmentLengths)
