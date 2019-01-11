import numpy as np
import sys
import helperFunctions
import copy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
try:
    import mpl_toolkits.mplot3d.axes3d as p3
except ImportError:
    print()
    pass


class chromophore:
    def __init__(self, chromoID, chromophoreCGSites, CGMorphologyDict, AAMorphologyDict, CGToAAIDMaster, parameterDict, simDims):
        self.ID = chromoID
        self.orcaInput = '/chromophores/inputORCA/single/%04d.inp' % (self.ID)
        self.orcaOutput = '/chromophores/outputORCA/single/%04d.out' % (self.ID)
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
        # The position of the chromophore can be calculated easily. Note that here, the `self.image' is the periodic image that the unwrapped_position of the chromophore is located in, relative to the original simulation volume.
        self.unwrappedPosn, self.posn, self.image = self.obtainChromophoreCoM(electronicallyActiveCGSites, flattenedCGToAAIDMaster, AAMorphologyDict, simDims)
        # A list of the important bonds for this chromophore from the morphology would be useful when determining
        # if a terminating group is already present on this monomer
        self.bonds = self.getImportantBonds(AAMorphologyDict['bond'])
        # Now to create a load of placeholder parameters to update later when we have the full list/energy levels
        # The self.neighbours list contains one element for each chromophore within parameterDict['maximumHopDistance']
        # of this one (including periodic boundary conditions). Its format is [[neighbour1ID, relativeImageOfNeighbour1],...]
        self.neighbours = []
        # Also find out the chromoIDs that are bonded directly to this one
        self.bondedChromos = []
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
        # The following lists are treated later on in the pipeline to add in additional neighbours to help model the
        # variable chromophore length in polymers, which permits the carrier to hop between pairs of monomers that are nearby,
        # but not necessarily the current monomer. Therefore, completeNeighbours shows the FINAL location of the charge 
        # carrier (no matter where it started). The completeNeighboursDeltaE and completeNeighboursTI are determined from the
        # conventional neighbour hopping for the nearby pairs. 
        self.completeNeighbours = []
        self.completeNeighboursDeltaE = []
        self.completeNeighboursTI = []
        pass

    def getImportantBonds(self, bondList):
        importantBonds = []
        for bond in bondList:
            if (bond[1] in self.AAIDs) and (bond[2] in self.AAIDs):
                importantBonds.append(bond)
        return importantBonds

    def obtainChromophoreCoM(self, electronicallyActiveCGSites, flattenedCGToAAIDMaster, AAMorphologyDict, simDims):
        # By using electronicallyActiveCGSites, determine the AAIDs for the electrically active proportion of the
        # chromophore, so that we can calculate its proper position. Again each element corresponds to each CG site
        # so the list needs to be flattened afterwards.
        electronicallyActiveAAIDs = [AAID for AAIDs in [flattenedCGToAAIDMaster[CGID] for CGID in electronicallyActiveCGSites] for AAID in AAIDs]
        electronicallyActiveUnwrappedPosns = [AAMorphologyDict['unwrapped_position'][AAID] for AAID in electronicallyActiveAAIDs]
        electronicallyActiveTypes = [AAMorphologyDict['type'][AAID] for AAID in electronicallyActiveAAIDs]
        # Calculate the chromophore's position in the morphology (CoM of all atoms in self.AAIDs from AAMorphologyDict)
        chromoUnwrappedPosn = helperFunctions.calcCOM(electronicallyActiveUnwrappedPosns, listOfAtomTypes=electronicallyActiveTypes)
        chromoWrappedPosn = copy.deepcopy(chromoUnwrappedPosn)
        chromoWrappedImage = [0, 0, 0]
        # Now calculate the wrapped position of the chromophore and its image
        for axis in range(3):
            simExtent = simDims[axis][1] - simDims[axis][0]
            while chromoWrappedPosn[axis] < simDims[axis][0]:
                chromoWrappedPosn[axis] += simExtent
                chromoWrappedImage[axis] -= 1
            while chromoWrappedPosn[axis] > simDims[axis][1]:
                chromoWrappedPosn[axis] -= simExtent
                chromoWrappedImage[axis] += 1
        return chromoUnwrappedPosn, chromoWrappedPosn, chromoWrappedImage

    def obtainElectronicSpecies(self, chromophoreCGSites, CGSiteTypes, CGToSpecies):
        electronicallyActiveSites = []
        currentChromophoreSpecies = None
        for CGSiteID in chromophoreCGSites:
            siteType = CGSiteTypes[CGSiteID]
            siteSpecies = CGToSpecies[siteType]
            if (siteSpecies != 'None'):
                if (currentChromophoreSpecies is not None) and (currentChromophoreSpecies != siteSpecies):
                    raise SystemError("PROBLEM - Multiple electronic species defined in the same chromophore. Please modify the chromophore generation code to fix this issue for your molecule!")
                else:
                    currentChromophoreSpecies = siteSpecies
                    electronicallyActiveSites.append(CGSiteID)
        return electronicallyActiveSites, currentChromophoreSpecies


def calculateChromophores(CGMorphologyDict, AAMorphologyDict, CGToAAIDMaster, parameterDict, simDims):
    # We make the assumption that a chromophore consists of one of each of the CG site types
    # described by the same template file. For instance, if we have 3 sites 'A', 'B' and 'C'
    # described in one file and one site 'D' described in another file then there are two
    # chromophores species described by A-B-C and D. This will be treated automatically
    # because the D's shouldn't be bonded to anything in the CGMorphologyDict if they are
    # small molecules.
    # Therefore, we need to assign each CG site in the morphology to a particular chromophore,
    # so first, it's important to generate a `neighbourlist' of all bonded atoms
    print()
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
    print()
    # Now let's create a list of all the chromophore instances which contain all of the
    # information we could ever want about them.
    chromophoreInstances = []
    for chromoID, chromophoreCGSites in chromophoreData.iteritems():
        print()
        sys.stdout.flush()
        chromophoreInstances.append(chromophore(chromoID, chromophoreCGSites, CGMorphologyDict, AAMorphologyDict, CGToAAIDMaster, parameterDict, simDims))
    print()
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


def determineNeighbours(chromophoreList, parameterDict, simDims, CGMorphologyDict):
    for chromophore1 in chromophoreList:
        print()
        sys.stdout.flush()
        for chromophore2 in chromophoreList:
            if chromophore1.ID == chromophore2.ID:
                continue
            deltaPosn = chromophore2.posn - chromophore1.posn
            relativeImageOfChromo2 = [0, 0, 0]
            # Consider periodic boundary conditions
            for axis in range(3):
                halfBoxLength = (simDims[axis][1] - simDims[axis][0]) / 2.0
                while deltaPosn[axis] > halfBoxLength:
                    deltaPosn[axis] -= simDims[axis][1] - simDims[axis][0]
                    relativeImageOfChromo2[axis] -= 1
                while deltaPosn[axis] < - halfBoxLength:
                    deltaPosn[axis] += simDims[axis][1] - simDims[axis][0]
                    relativeImageOfChromo2[axis] += 1
            separation = np.linalg.norm(deltaPosn)
            # If proximity is within tolerance, add these chromophores as neighbours
            if separation <= parameterDict['maximumHopDistance']:
                # Only add the neighbours if they haven't already been added so far
                chromo1NeighbourIDs = [neighbourData[0] for neighbourData in chromophore1.neighbours]
                chromo2NeighbourIDs = [neighbourData[0] for neighbourData in chromophore2.neighbours]
                # Determine if these chromophores are bonded
                for bond in CGMorphologyDict['bond']:
                    if ((bond[1] in chromophore1.CGIDs) and (bond[2] in chromophore2.CGIDs)) or ((bond[2] in chromophore1.CGIDs) and (bond[1] in chromophore2.CGIDs)):
                        if chromophore1.ID not in chromophore2.bondedChromos:
                            chromophore2.bondedChromos.append(chromophore1.ID)
                        if chromophore2.ID not in chromophore1.bondedChromos:
                            chromophore1.bondedChromos.append(chromophore2.ID)
                        break
                # Also, make the deltaE and the Tij lists as long as the neighbour lists for easy access later
                if chromophore2.ID not in chromo1NeighbourIDs:
                    chromophore1.neighbours.append([chromophore2.ID, relativeImageOfChromo2])
                    chromophore1.neighboursDeltaE.append(None)
                    chromophore1.neighboursTI.append(None)
                if chromophore1.ID not in chromo2NeighbourIDs:
                    chromophore2.neighbours.append([chromophore1.ID, list(-np.array(relativeImageOfChromo2))])
                    chromophore2.neighboursDeltaE.append(None)
                    chromophore2.neighboursTI.append(None)
        # DEBUG TESTING
        # if chromophore1.ID == 1961:
        #     print ""
        #     print chromophore1.posn
        #     print chromophore1.AAIDs
        #     for neighbour in chromophore1.neighbours:
        #         print neighbour[0], neighbour[1], chromophoreList[neighbour[0]].posn
        #     plotChromoNeighbours(chromophore1, chromophoreList, simDims)
    print()
    return chromophoreList


def plotChromoNeighbours(chromophore, chromophoreList, simDims):
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    ax.scatter(chromophore.posn[0], chromophore.posn[1], chromophore.posn[2], s=50, c='r')
#    for chromophore in chromophoreList:
#        position = chromophore.posn
#        for axis in range(3):
#            while position[axis] >= simDims[axis][1]:
#                position[axis] -= simDims[axis][1] - simDims[axis][0]
#            while position[axis] <= simDims[axis][0]:
#                position[axis] += simDims[axis][1] - simDims[axis][0]
#        ax.scatter(position[0], position[1], position[2], s = 5, c = 'k')
    print()
    print()
    for neighbour in chromophore.neighbours:
        neighbourID = neighbour[0]
        neighbourImage = neighbour[1]
        neighbourChromo = chromophoreList[neighbourID]
        if neighbourChromo.ID != neighbourID:
            raise SystemError("WRONG CHROMO")
        neighbourPosn = copy.deepcopy(neighbourChromo.posn)
        for axis in range(3):
            simLength = simDims[axis][1] - simDims[axis][0]
            neighbourPosn[axis] += neighbourImage[axis] * simLength
        print()
        ax.scatter(neighbourPosn[0], neighbourPosn[1], neighbourPosn[2], s=50, c='b')
        ax.scatter(neighbourChromo.posn[0], neighbourChromo.posn[1], neighbourChromo.posn[2], s=50, c='g')
    # Finally, draw box
    ax.set_xlim(1.1 * np.array(simDims[0]))
    ax.set_ylim(1.1 * np.array(simDims[1]))
    ax.set_zlim(1.1 * np.array(simDims[2]))
    plt.savefig('./test.pdf')
    plt.show()
    exit()


def execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList):
    simDims = [[-AAMorphologyDict['lx'] / 2.0, AAMorphologyDict['lx'] / 2.0], [-AAMorphologyDict['ly'] / 2.0, AAMorphologyDict['ly'] / 2.0], [-AAMorphologyDict['lz'] / 2.0, AAMorphologyDict['lz'] / 2.0]]
    chromophoreList = calculateChromophores(CGMorphologyDict, AAMorphologyDict, CGToAAIDMaster, parameterDict, simDims)
    chromophoreList = determineNeighbours(chromophoreList, parameterDict, simDims, CGMorphologyDict)
    # Now we have updated the chromophoreList, rewrite the pickle with this new information.
    pickleName = parameterDict['outputDir'] + '/' + parameterDict['morphology'][:-4] + '/code/' + parameterDict['morphology'][:-4] + '.pickle'
    helperFunctions.writePickle((AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList), pickleName)
    return AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList


if __name__ == "__main__":
    try:
        pickleFile = sys.argv[1]
    except:
        print()
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList = helperFunctions.loadPickle(pickleFile)
    execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList)
