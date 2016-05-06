import numpy as np
import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cPickle as pickle
import cme_utils
import helperFunctions
import chromophores
try:
    import mpl_toolkits.mplot3d.axes3d as p3
except ImportError:
    print "Could not import 3D plotting engine, calling the plotMolecule3D function will result in an error!"
    pass


def getFunctionalGroups(molecule, CGtoAAIDs, CGBonds):
    monomerData = []
    for thioID in sorted(CGtoAAIDs.keys()):
        if CGtoAAIDs[thioID][0] == 'thio':
            currentMonomer = [CGtoAAIDs[thioID][1]]
            for bond in CGBonds:
                if (bond[0] == 'bondB') and (bond[1] == thioID):
                    alk1ID = bond[2]
                    currentMonomer.append(CGtoAAIDs[alk1ID][1])
                    break
                elif (bond[0] == 'bondB') and (bond[2] == thioID):
                    alk1ID = bond[1]
                    currentMonomer.append(CGtoAAIDs[alk1ID][1])
                    break
            for bond in CGBonds:
                if (bond[0] == 'bondC') and (bond[1] == alk1ID):
                    alk2ID = bond[2]
                    currentMonomer.append(CGtoAAIDs[alk2ID][1])
                    break
                elif (bond[0] == 'bondC') and (bond[2] == alk1ID):
                    alk2ID = bond[1]
                    currentMonomer.append(CGtoAAIDs[alk2ID][1])
                    break
            monomerData.append(currentMonomer)
        else:
            continue
    moleculeEnds = []
    for monomerNo, monomer in enumerate(monomerData):
        if len(monomer[0]) == 7:
            moleculeEnds.append(monomerNo)
    monomerData = np.array(monomerData)
    thioRings = monomerData[:,0]
    alk1Groups = monomerData[:,1]
    alk2Groups = monomerData[:,2]
    return thioRings, alk1Groups, alk2Groups, moleculeEnds


def obtainBackboneData(molecule, thioRings):
    thioGroups = []
    moleculeEnds = []
    for thioRing in thioRings:
        thioGroup = {'index': [], 'type': [], 'mass': [], 'position': [], 'unwrapped_position': []}
        # To be of the form, [IDs], [Types], [Masses], [Positions]
        # The S and the preceeding 4 atoms make up the thiophene ring,
        # along with the following H (2 Hs if at the end of the molecule)
        for atomID in thioRing:
            thioGroup['index'].append(atomID)
            thioGroup['type'].append(molecule['type'][atomID])
            thioGroup['mass'].append(molecule['mass'][atomID])
            thioGroup['unwrapped_position'].append(molecule['unwrapped_position'][atomID])
            thioGroup['position'].append(molecule['position'][atomID])
        thioGroup['COM'] = helperFunctions.calcCOM(thioGroup['position'], thioGroup['mass'])
        thioGroup['unwrapped_COM'] = helperFunctions.calcCOM(thioGroup['unwrapped_position'], thioGroup['mass'])
        # Plane of the thiophene ring is given by the vector between the
        # C1 and C10 (atoms 0 and 3 in the ring)
        thioPlaneVector = helperFunctions.findAxis(molecule['unwrapped_position'][thioRing[0]], molecule['unwrapped_position'][thioRing[3]])
        thioGroup['plane'] = thioPlaneVector
        # Normal to plane of thiophene ring is given by the cross product
        # between the C1-C10 vector and the C1-C2 vector
        thioNormalVector = helperFunctions.normaliseVec(np.cross(thioPlaneVector, helperFunctions.findAxis(molecule['unwrapped_position'][thioRing[0]], molecule['unwrapped_position'][thioRing[1]])))
        thioGroup['normal'] = thioNormalVector
        thioGroups.append(thioGroup)
    return thioGroups


def calculatePersistenceLength(molecule, molName, outputDir, moleculeBackbone, moleculeEnds):
    distances = []
    cosThetas = []
    for monomerID1, thioRing1 in enumerate(moleculeBackbone):
        for monomerID2, thioRing2 in enumerate(moleculeBackbone):
            #separation = helperFunctions.calculateSeparation(thioRing1['COM'], thioRing2['COM'])
            # Separation in monomer units can be given by monomerID2 - monomerID1
            distances.append(monomerID2 - monomerID1)
            cosThetas.append(np.dot(thioRing1['plane'], thioRing2['plane']))
    #distances, cosThetas = helperFunctions.parallelSort(distances, cosThetas)
    cosThetaDictionary = {}
    for index, monomerDistance in enumerate(distances):
        if abs(monomerDistance) not in cosThetaDictionary:
            cosThetaDictionary[abs(monomerDistance)] = []
        cosThetaDictionary[abs(monomerDistance)].append(cosThetas[index])


    distanceToPlot = []
    cosThetaToPlot = []
    ybars = []
        
    for monomerDistance in cosThetaDictionary:
        distanceToPlot.append(monomerDistance)
        cosThetaToPlot.append(np.average(cosThetaDictionary[monomerDistance]))
        ybars.append(np.std(cosThetaDictionary[monomerDistance])/np.sqrt(len(cosThetaDictionary[monomerDistance])))
    autoCorrelationArray = cme_utils.analyze.autocorr.autocorr1D(cosThetaToPlot)
    persistenceLength = helperFunctions.linearInterpDescendingY(0, distanceToPlot, autoCorrelationArray)
    if persistenceLength == None:
        print "All lengths within this molecule are correlated. Setting persistence length to maximum ("+str(np.max(distanceToPlot[:len(autoCorrelationArray)]))+")..."
        persistenceLength = np.max(distanceToPlot[:len(autoCorrelationArray)])
    #else:
        #print "Persistence Length for this molecule =", persistenceLength, "monomers."
    plt.figure()
    plt.plot(distanceToPlot, cosThetaToPlot, 'ro')
    plt.errorbar(distanceToPlot, cosThetaToPlot, yerr=ybars, linestyle='None')
    plt.title("PL = "+str(persistenceLength))
    plt.savefig(outputDir+'/molecules/'+molName.replace('.POSCAR', '_PL.png'))
    plt.close()
    return persistenceLength


def plotMolecule3D(molecule, moleculeBackbone):
    allThioAtoms = []
    for thioRing in moleculeBackbone:
        allThioAtoms += list(np.array(thioRing['position']))
    allThioAtoms = np.array(allThioAtoms)
    allAtoms = np.array(molecule['unwrapped_position'])
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
    ax.scatter(allThioAtoms[:,0], allThioAtoms[:,1], allThioAtoms[:,2], s = 20, c = 'k')
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
    bendingAngleTolerance = np.pi#0.0#np.pi#*(1/3.)
    # The bending angle, beyond which the pi conjugation is said to be broken
    torsionAngleTolerance = np.pi#0.0#np.pi#*(3/8.) # The torsional angle, beyond which, the pi conjugation is said to be broken
    previousThioPlaneAxis = moleculeBackbone[0]['plane']
    previousThioNormalAxis = moleculeBackbone[0]['normal']
    bendingAngles = []
    torsionAngles = []
    chromophoreIDs = [[0]]
    thioCOMs = [[moleculeBackbone[0]['COM']]]
    unwrappedThioCOMs = [[moleculeBackbone[0]['unwrapped_COM']]]
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
            #print "Torsion has broken conjugation between", thioNumber-1, "and", thioNumber, ":", torsionAngle*180/np.pi, ">", torsionAngleTolerance*180/np.pi
            torsionOverThreshold += 1
            conjugationBroken = True
        if (bendingAngle > bendingAngleTolerance):
            #print "Bending has broken conjugation between", thioNumber-1, "and", thioNumber, ":", bendingAngle*180/np.pi, ">", bendingAngleTolerance*180/np.pi
            bendingOverThreshold += 1
            conjugationBroken = True
        if conjugationBroken == True:
            chromophoreIDs.append([])
            thioCOMs.append([])
            unwrappedThioCOMs.append([])
            conjugationBroken = False
        chromophoreIDs[-1].append(thioNumber)
        thioCOMs[-1].append(thioRing['COM'])
        unwrappedThioCOMs[-1].append(thioRing['unwrapped_COM'])
        bendingAngles.append(bendingAngle)
        torsionAngles.append(torsionAngle)
        previousThioPlaneAxis = np.copy(currentThioPlaneAxis)
        previousThioNormalAxis = np.copy(currentThioNormalAxis)
    # print "\nMolecule", molName
    # print "The torsional angle was over the threshold ("+str(torsionAngleTolerance)+")", torsionOverThreshold, "times."
    # print "The bending angle was over the threshold ("+str(bendingAngleTolerance)+")", bendingOverThreshold, "times."
    #plotHist(bendingAngles, outputDir+'/molecules/'+molName.replace('.POSCAR', '_Bend.png'), angle=True)
    #plotHist(torsionAngles, outputDir+'/molecules/'+molName.replace('.POSCAR', '_Tor.png'), angle=True)
    return bendingAngles, torsionAngles, chromophoreIDs, thioCOMs, unwrappedThioCOMs


def plotHist(data, outputFile, bins=20, angle=False):
    if angle == True:
        data=np.array(data)*180/np.pi
        bins = np.arange(0,90,5)
    plt.figure()
    plt.hist(data, bins, normed=1)
    plt.savefig(outputFile)
    plt.close()

def execute(morphologyFile, AAfileName, CGMoleculeDict, UnrelaxedAAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize):
    morphologyName = morphologyFile[helperFunctions.findIndex(morphologyFile,'/')[-1]+1:]
    outputDir = './outputFiles'
    morphologyList = os.listdir(outputDir)
    for allMorphologies in morphologyList:
        if morphologyName in allMorphologies:
            outputDir += '/'+morphologyName
            break
    CGBonds = []
    for CGBond in CGMoleculeDict['bond']:
        if CGBond not in CGBonds:
            CGBonds.append(CGBond)

    
    # The current AAMorphologyDict has been read in from the Pickle File which means it's pre-phase 1.
    # Need to update it to the latest by reading in the relaxed XML
    slashList = helperFunctions.findIndex(AAfileName, '/')
    inputFileName = AAfileName[:slashList[-1]+1]+'relaxed_'+AAfileName[slashList[-1]+1:]
    print "Loading relaxed morphology data..."
    AAMorphologyDict = helperFunctions.loadMorphologyXML(inputFileName)
    inverseSScale = helperFunctions.getsScale(outputDir, morphologyName)
    # Scale up the morphology (needs to be in Angstroems for ORCA so do it now)
    AAMorphologyDict = helperFunctions.scale(AAMorphologyDict, inverseSScale)
    AAMorphologyDict = helperFunctions.addUnwrappedPositions(AAMorphologyDict)
    ## Work only on the unwrapped positions
    #AAMorphologyDict['position'] = AAMorphologyDict['unwrapped_position']
    

    
    # for molID, molAAIDs in enumerate(moleculeAAIDs):
    #     moleculeDict = helperFunctions.loadDict(AAMorphologyDict, molAAIDs)
    #     moleculeList.append(moleculeDict)
    #     molName = str(molID)
    #     while len(molName) < 4:
    #         molName = '0'+molName
    #     molNames.append('mol'+molName+'.POSCAR')
    # # IGNORE: CANNOT DO THIS, LOADING FROM POSCAR SCREWS UP POSITIONING
    # for poscarFile in poscarList:
    #     if '.POSCAR' in poscarFile:
    #         moleculeDict = helperFunctions.loadPoscar(outputDir+'/molecules/'+poscarFile)
    #         moleculeDict = helperFunctions.addMasses(moleculeDict)
    #         moleculeList.append(moleculeDict)
    #         molNames.append(poscarFile)

    moleculeAAIDs = helperFunctions.getAAIDsByMolecule(CGtoAAIDs)
    moleculeList = []
    molNames = []
    print "Obtaining molecule dictionaries..."
    for molNo, molecule in enumerate(moleculeAAIDs):
        print "Loading molecule", molNo+1, "of", str(len(moleculeAAIDs))+"...\r",
        molName = str(molNo)
        while len(molName) < 3:
            molName = '0'+molName
        molNames.append('mol'+molName+'.POSCAR')
        moleculeDict = helperFunctions.loadDict(AAMorphologyDict, molecule, outputDir+'/molecules/mol'+molName+'.pickle')
        moleculeList.append(moleculeDict)
    print "\n"
    
    morphologyBackboneData = [] # Dictionaries of molecule data for each molecule in the system
    rollingAtomID = 0
    for molNo, molecule in enumerate(moleculeList):
        # if molNo != 69:
        #     continue
        print "Examining molecule number", molNo+1, "of", str(len(moleculeList))+"...\r",
        thioRings, alk1Groups, alk2Groups, moleculeEnds = getFunctionalGroups(molecule, CGtoAAIDs[0], CGBonds)
        moleculeBackbone = obtainBackboneData(molecule, thioRings)
        endToEndDistance = helperFunctions.calculateSeparation(moleculeBackbone[moleculeEnds[0]]['COM'], moleculeBackbone[moleculeEnds[1]]['COM'])
        persistenceLength = calculatePersistenceLength(molecule, molNames[molNo], outputDir, moleculeBackbone, moleculeEnds)
        bendingAngles, torsionAngles, chromophoreIDs, thiopheneCOMs, unwrappedThiopheneCOMs = calculateChromophores(molecule, molNames[molNo], outputDir, moleculeBackbone)
        morphologyChromophores = []
        thioCOMs = []
        for chromophore in chromophoreIDs:
            morphologyChromophores.append([])
            for monomerID in chromophore:
                morphologyChromophores[-1] += list(np.array(thioRings[monomerID])+rollingAtomID)+list(np.array(alk1Groups[monomerID])+rollingAtomID)+list(np.array(alk2Groups[monomerID])+rollingAtomID)
        moleculeDictionary = {'ends': moleculeEnds, 'endToEndDistance': endToEndDistance, 'persistenceLength': persistenceLength, 'bendingAngles': bendingAngles, 'torsionAngles': torsionAngles, 'chromophores': chromophoreIDs, 'morphologyChromophores': morphologyChromophores, 'thioCOMs': thiopheneCOMs, 'unwrappedThioCOMs': unwrappedThiopheneCOMs}
        morphologyBackboneData.append(moleculeDictionary)
        rollingAtomID += len(molecule['type'])
        # plotMolecule3D(molecule, moleculeBackbone)
    print "\n"
    print "Plotting morphology histograms..."
    endToEndDistance = []
    bendingAngles = []
    torsionAngles = []
    persistenceLengths = []
    for molDict in morphologyBackboneData:
        endToEndDistance.append(molDict['endToEndDistance'])
        bendingAngles += molDict['bendingAngles']
        torsionAngles += molDict['torsionAngles']
        persistenceLengths.append(molDict['persistenceLength'])
    plotHist(endToEndDistance, outputDir+'/morphology/EndToEndDistances.png')
    plotHist(bendingAngles, outputDir+'/morphology/BendingAngles.png', angle=True)
    plotHist(torsionAngles, outputDir+'/morphology/TorsionAngles.png', angle=True)
    plotHist(persistenceLengths, outputDir+'/morphology/PersistenceLengths.png', angle=False)
    segmentLengths = []
    for molDict in morphologyBackboneData:
        for segment in molDict['chromophores']:
            segmentLengths.append(len(segment))
    print "\n--== MORPHOLOGY STATISTICS ==--"
    print "Average end-to-end distance in angstroms (this has been scaled correctly) =", np.average(endToEndDistance)
    print "Average bending angle in degrees =", np.average(bendingAngles)*180/np.pi
    print "Average torsional angle in degrees =", np.average(torsionAngles)*180/np.pi
    print "Average persistence length in monomer units =", np.average(persistenceLengths)
    print "Average chromophore length in monomers =", np.average(segmentLengths), "\n"

    # NOW, WRITE THE CHROMOPHORES OUT TO AN ORCA INPUT FILE OR SOMETHING IN ./outputFiles/<morphName>/chromophores/single.
    chromophores.obtain(AAMorphologyDict, morphologyBackboneData, boxSize, outputDir, moleculeAAIDs)
    # Also, get their COM posn so that we can split them into neighbouring pairs and store it in ../pairs
    return morphologyFile, AAfileName, CGMoleculeDict, UnrelaxedAAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize
    
    
def loadPickle(morphologyFile):
    morphologyName = morphologyFile[helperFunctions.findIndex(morphologyFile,'/')[-1]+1:]
    outputDir = './outputFiles'
    morphologyList = os.listdir(outputDir)
    for allMorphologies in morphologyList:
        if morphologyName in allMorphologies:
            outputDir += '/'+morphologyName
            break
    pickleFound = False
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
        (AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize) = pickle.load(pickleFile)
    morphologyFile, AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize = execute(morphologyFile, AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize)
    return morphologyFile, AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize


if __name__ == '__main__':
    morphologyFile = sys.argv[1]
    loadPickle(morphologyFile)
