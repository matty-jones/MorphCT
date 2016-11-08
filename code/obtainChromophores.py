import numpy as np
import sys
import os
#import matplotlib
#matplotlib.use('Agg')
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




class chromophore:
    def __init__(self, chromoID, chromophoreCGSites, CGMorphologyDict, AAMorphologyDict, CGToAAIDMaster, parameterDict):
        self.ID = chromoID
        self.orcaInput = parameterDict['outputDir']+'/'+parameterDict['morphology'][:-4]+'/chromophores/inputORCA/single/%04d.inp' % (self.ID)
        self.orcaOutput = parameterDict['outputDir']+'/'+parameterDict['morphology'][:-4]+'/chromophores/outputORCA/single/%04d.out' % (self.ID)
        self.CGIDs = chromophoreCGSites
        # Determine whether this chromophore is a donor or an acceptor, as well as the site types that have been
        # defined as the electronically active in the chromophore
        electronicallyActiveCGSites, self.species = self.obtainElectronicSpecies(chromophoreCGSites, CGMorphologyDict['type'], parameterDict['CGSiteSpecies'])
        # CGToAAIDMaster is a list of dictionaries where each list element corresponds to a new molecule.
        # Firstly, flatten this out so that it becomes a single CG:AAID dictionary
        flattenedCGToAAIDMaster = {dictKey: dictVal[1] for dictionary in CGToAAIDMaster for dictKey, dictVal in dictionary.iteritems()}
        # Now, using chromophoreCGSites as the keys, build up a list of all of the AAIDs in the chromophore,
        # where each element corresponds to each CG site, and then flatten it.
        self.AAIDs = [AAID for AAIDs in [flattenedCGToAAIDMaster[CGID] for CGID in chromophoreCGSites] for AAID in AAIDs]
        # The position of the chromophore can be calculated easily
        self.posn = self.obtainChromophoreCoM(electronicallyActiveCGSites, flattenedCGToAAIDMaster, AAMorphologyDict)
        # A list of the important bonds for this chromophore from the morphology would be useful when determining
        # if a terminating group is already present on this monomer
        self.bonds = self.getImportantBonds(AAMorphologyDict['bond'])
        # Now to create a load of placeholder parameters to update later when we have the full list/energy levels
        # The self.neighbours list contains one element for each chromophore within parameterDict['maximumHopDistance']
        # of this one (including periodic boundary conditions). Its format is [[neighbour1ID, relativeImageOfNeighbour1],...]
        self.neighbours = []
        # The molecular orbitals of this chromophore have not yet been calculated, but they will simply be floats.
        self.HOMO = None
        self.HOMO_1 = None
        self.LUMO = None
        self.LUMO_1 = None
        # The neighbourDeltaE and neighbourTI are lists where each element describes the different in important molecular
        # orbital or transfer integral between this chromophore and each neighbour. The list indices here are the same as
        # in self.neighbours for coherence.
        self.neighboursDeltaE = []
        self.neighboursTI = []
        # Placeholder parameters to be updated in the pickle later?:
            # HOMO, HOMO+1, LUMO, LUMO-1
            # Transfer Integral List to each neighbour
            # Number of times particular pathways have been hopped along? So like chromoID1-chromoID2 = 27 which increments each time a carrier hops from chromoID1 to chromoID2 
        pass

    def getImportantBonds(self, bondList):
        importantBonds = []
        for bond in bondList:
            if (bond[1] in self.AAIDs) and (bond[2] in self.AAIDs):
                importantBonds.append(bond)
        return importantBonds

    def obtainChromophoreCoM(self, electronicallyActiveCGSites, flattenedCGToAAIDMaster, AAMorphologyDict):
        # By using electronicallyActiveCGSites, determine the AAIDs for the electrically active proportion of the
        # chromophore, so that we can calculate its proper position. Again each element corresponds to each CG site
        # so the list needs to be flattened afterwards.
        electronicallyActiveAAIDs = [AAID for AAIDs in [flattenedCGToAAIDMaster[CGID] for CGID in electronicallyActiveCGSites] for AAID in AAIDs]
        electronicallyActivePosns = [AAMorphologyDict['unwrapped_position'][AAID] for AAID in electronicallyActiveAAIDs]
        electronicallyActiveTypes = [AAMorphologyDict['type'][AAID] for AAID in electronicallyActiveAAIDs]
        # Calculate the chromophore's position in the morphology (CoM of all atoms in self.AAIDs from AAMorphologyDict)
        return helperFunctions.calcCOM(electronicallyActivePosns, listOfAtomTypes = electronicallyActiveTypes)

    def obtainElectronicSpecies(self, chromophoreCGSites, CGSiteTypes, CGToSpecies):
        electronicallyActiveSites = []
        currentChromophoreSpecies = None
        for CGSiteID in chromophoreCGSites:
            siteType = CGSiteTypes[CGSiteID]
            siteSpecies = CGToSpecies[siteType]
            if (siteSpecies != 'None'):
                if (currentChromophoreSpecies != None) and (currentChromophoreSpecies != siteSpecies):
                    raise SystemError("PROBLEM - Multiple electronic species defined in the same chromophore. Please modify the chromophore generation code to fix this issue for your molecule!")
                else:
                    currentChromophoreSpecies = siteSpecies
                    electronicallyActiveSites.append(CGSiteID)
        return electronicallyActiveSites, currentChromophoreSpecies


def calculateChromophores(CGMorphologyDict, AAMorphologyDict, CGToAAIDMaster, parameterDict):
    # We make the assumption that a chromophore consists of one of each of the CG site types
    # described by the same template file. For instance, if we have 3 sites 'A', 'B' and 'C'
    # described in one file and one site 'D' described in another file then there are two 
    # chromophores species described by A-B-C and D. This will be treated automatically 
    # because the D's shouldn't be bonded to anything in the CGMorphologyDict if they are 
    # small molecules.
    # Therefore, we need to assign each CG site in the morphology to a particular chromophore,
    # so first, it's important to generate a `neighbourlist' of all bonded atoms
    print "Determining chromophores in the system..."
    bondedAtoms = helperFunctions.obtainBondedList(CGMorphologyDict['bond'])
    chromophoreList = [i for i in range(len(CGMorphologyDict['type']))]
    for CGSiteID, chromophoreID in enumerate(chromophoreList):
        CGSiteType = CGMorphologyDict['type'][CGSiteID]
        typesInThisChromophore = [CGSiteType]
        chromophoreList, typesInThisChromophore = updateChromophores(CGSiteID, chromophoreList, bondedAtoms, CGMorphologyDict['type'], typesInThisChromophore, parameterDict)
    chromophoreData = {}
    for atomID, chromoID in enumerate(chromophoreList):
        if chromoID not in chromophoreData.keys():
            chromophoreData[chromoID] = [atomID]
        else:
            chromophoreData[chromoID].append(atomID)
    # Now rename the chromophore IDs so that they increment sensibly (they will be used later 
    # for the ORCA files)
    oldKeys = sorted(chromophoreData.keys())
    for newKey, oldKey in enumerate(oldKeys):
        chromophoreData[newKey] = chromophoreData.pop(oldKey)
    print str(len(chromophoreData.keys()))+" chromophores successfully identified!"
    # Now let's create a list of all the chromophore instances which contain all of the 
    # information we could ever want about them.
    chromophoreInstances = []
    for chromoID, chromophoreCGSites in chromophoreData.iteritems():
        print "\rCalculating properties of chromophore %04d of %04d..." % (chromoID, len(chromophoreData.keys())-1),
        sys.stdout.flush()
        chromophoreInstances.append(chromophore(chromoID, chromophoreCGSites, CGMorphologyDict, AAMorphologyDict, CGToAAIDMaster, parameterDict))
    print ""
    return chromophoreInstances


def updateChromophores(atomID, chromophoreList, bondedAtoms, CGTypeList, typesInThisChromophore, parameterDict):
    # Recursively add all neighbours of atom number atomID to this chromophore, providing the same type does not already exist in it
    try:
        for bondedAtom in bondedAtoms[atomID]:
            bondedType = CGTypeList[bondedAtom]
            # First check that the bondedAtom's type is not already in this chromophore.
            # Also, check that the type to be added is of the same electronic species as the ones added previously, or == 'None'
            if (bondedType not in typesInThisChromophore) and ((parameterDict['CGSiteSpecies'][bondedType] == 'None') or (parameterDict['CGSiteSpecies'][bondedType] == list(set([parameterDict['CGSiteSpecies'][x] for x in typesInThisChromophore]))[0])):
                # If the atomID of the bonded atom is larger than that of the current one,
                # update the bonded atom's ID to the current one's to put it in this chromophore,
                # then iterate through all of the bonded atom's neighbours
                if chromophoreList[bondedAtom] > chromophoreList[atomID]:
                    chromophoreList[bondedAtom] = chromophoreList[atomID]
                    typesInThisChromophore.append(bondedType)
                    chromophoreList, typesInThisChromophore = updateChromophores(bondedAtom, chromophoreList, bondedAtoms, CGTypeList, typesInThisChromophore, parameterDict)
                # If the atomID of the current atom is larger than that of the bonded one,
                # update the current atom's ID to the bonded one's to put it in this chromophore,
                # then iterate through all of the current atom's neighbours
                elif chromophoreList[bondedAtom] < chromophoreList[atomID]:
                    chromophoreList[atomID] = chromophoreList[bondedAtom]
                    typesInThisChromophore.append(CGTypeList[atomID])
                    chromophoreList, typesInThisChromophore = updateChromophores(atomID, chromophoreList, bondedAtoms, CGTypeList, typesInThisChromophore, parameterDict)
                # Else: both the current and the bonded atom are already known to be in this 
                # chromophore, so we don't have to do anything else.
    except KeyError:
        # This means that there are no bonded CG sites (i.e. it's a single chromophore)
        pass
    return chromophoreList, typesInThisChromophore


def determineNeighbours(chromophoreList, parameterDict, simDims):
    for chromophore1 in chromophoreList:
        print "\rIdentifying neighbours of chromophore %04d of %04d..." % (chromophore1.ID, len(chromophoreList)-1),
        sys.stdout.flush()
        for chromophore2 in chromophoreList:
            if chromophore1.ID == chromophore2.ID:
                continue
            deltaPosn = chromophore2.posn - chromophore1.posn
            relativeImageOfChromo2 = [0, 0, 0]
            # Consider periodic boundary conditions
            for axis in range(3):
                while deltaPosn[axis] > simDims[axis][1]:
                    deltaPosn[axis] -= simDims[axis][1]-simDims[axis][0]
                    relativeImageOfChromo2[axis] += 1
                while deltaPosn[axis] < simDims[axis][0]:
                    deltaPosn[axis] += simDims[axis][1]-simDims[axis][0]
                    relativeImageOfChromo2[axis] -= 1
            separation = np.linalg.norm(deltaPosn)
            # If proximity is within tolerance, add these chromophores as neighbours
            if separation <= parameterDict['maximumHopDistance']:
                chromophore1.neighbours.append([chromophore2.ID, relativeImageOfChromo2])
                chromophore2.neighbours.append([chromophore1.ID, list(-np.array(relativeImageOfChromo2))])
        # plotChromoNeighbours(chromophore1, chromophoreList, simDims)
    print ""
    return chromophoreList


def plotChromoNeighbours(chromophore, chromophoreList, simDims):
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    ax.scatter(chromophore.posn[0], chromophore.posn[1], chromophore.posn[2], s = 50, c = 'r')
    for neighbour in chromophore.neighbours:
        neighbourID = neighbour[0]
        neighbourImage = neighbour[1]
        neighbourChromo = chromophoreList[neighbourID]
        if neighbourChromo.ID != neighbourID:
            raise SystemError("WRONG CHROMO")
        neighbourPosn = neighbourChromo.posn
        for axis in range(3):
            simLength = simDims[axis][1] - simDims[axis][0]
            while neighbourImage[axis] > 0:
                neighbourImage[axis] -= 1
                neighbourPosn[axis] -= simLength
            while neighbourImage[axis] < 0:
                neighbourImage[axis] += 1
                neighbourPosn[axis] += simLength
        ax.scatter(neighbourPosn[0], neighbourPosn[1], neighbourPosn[2], s = 50, c = 'b')
    print "Showing..."
    plt.savefig('./test.pdf')


def execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList):
    chromophoreList = calculateChromophores(CGMorphologyDict, AAMorphologyDict, CGToAAIDMaster, parameterDict)
    chromophoresToPlot = range(15)
    #plotMolecule3D(chromophoresToPlot, chromophoreList, [[-AAMorphologyDict['lx']/2.0, AAMorphologyDict['lx']/2.0], [-AAMorphologyDict['ly']/2.0, AAMorphologyDict['ly']/2.0], [-AAMorphologyDict['lz']/2.0, AAMorphologyDict['lz']/2.0]])
    simDims = [[-AAMorphologyDict['lx']/2.0, AAMorphologyDict['lx']/2.0], [-AAMorphologyDict['ly']/2.0, AAMorphologyDict['ly']/2.0], [-AAMorphologyDict['lz']/2.0, AAMorphologyDict['lz']/2.0]]
    chromophoreList = determineNeighbours(chromophoreList, parameterDict, simDims)
    # Now we have updated the chromophoreList, rewrite the pickle with this new information.
    pickleName = parameterDict['outputDir']+'/'+parameterDict['morphology'][:-4]+'/code/'+parameterDict['morphology'][:-4]+'.pickle'
    helperFunctions.writePickle((AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList), pickleName)
    return AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList


if __name__ == "__main__":
    try:
        pickleFile = sys.argv[1]
    except:
        print "Please specify the pickle file to load to continue the pipeline from this point."
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList = helperFunctions.loadPickle(pickleFile)
    execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList)
