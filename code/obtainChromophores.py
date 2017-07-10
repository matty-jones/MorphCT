import numpy as np
import sys
import helperFunctions
import copy


class chromophore:
    def __init__(self, chromoID, chromophoreCGSites, CGMorphologyDict, AAMorphologyDict, CGToAAIDMaster, parameterDict, simDims):
        self.ID = chromoID
        self.orcaInput = '/chromophores/inputORCA/single/%04d.inp' % (self.ID)
        self.orcaOutput = '/chromophores/outputORCA/single/%04d.out' % (self.ID)
        self.CGIDs = chromophoreCGSites
        # Determine whether this chromophore is a donor or an acceptor, as well as the site types that have been
        # defined as the electronically active in the chromophore
        if CGMorphologyDict is not None:
            # Normal operation
            self.CGTypes = sorted(list(set([CGMorphologyDict['type'][CGID] for CGID in self.CGIDs])))
            electronicallyActiveCGSites, self.species = self.obtainElectronicSpecies(chromophoreCGSites, CGMorphologyDict['type'], parameterDict['CGSiteSpecies'])
            # CGToAAIDMaster is a list of dictionaries where each list element corresponds to a new molecule.
            # Firstly, flatten this out so that it becomes a single CG:AAID dictionary
            flattenedCGToAAIDMaster = {dictKey: dictVal[1] for dictionary in CGToAAIDMaster for dictKey, dictVal in dictionary.items()}
            # Now, using chromophoreCGSites as the keys, build up a list of all of the AAIDs in the chromophore,
            # where each element corresponds to each CG site, and then flatten it.
            self.AAIDs = [AAID for AAIDs in [flattenedCGToAAIDMaster[CGID] for CGID in chromophoreCGSites] for AAID in AAIDs]
            # By using electronicallyActiveCGSites, determine the AAIDs for the electrically active proportion of the
            # chromophore, so that we can calculate its proper position. Again each element corresponds to each CG site
            # so the list needs to be flattened afterwards.
            electronicallyActiveAAIDs = [AAID for AAIDs in [flattenedCGToAAIDMaster[CGID] for CGID in electronicallyActiveCGSites] for AAID in AAIDs]
        else:
            # No fine-graining has been performed by MorphCT, so we know that the input morphology is already atomistic.
            if len(parameterDict['CGSiteSpecies']) == 1:
                # If the morphology contains only a single type of electronic species, then the parameterDict['CGSiteSpecies'] should only have one entry, and we can set all chromophores to be this species.
                electronicallyActiveCGSites = chromophoreCGSites
                electronicallyActiveAAIDs = chromophoreCGSites
                self.species = list(parameterDict['CGSiteSpecies'].values())[0]
            elif (len(parameterDict['CGSiteSpecies']) == 0) and (len(parameterDict['AARigidBodySpecies']) > 0):
                # If the CGSiteSpecies have not been specified, then look to the AARigidBodySpecies dictionary to determine which rigid bodies are donors and which are acceptors
                electronicallyActiveAAIDs = []
                for AAID in chromophoreCGSites:
                    if AAMorphologyDict['body'][AAID] != -1:
                        electronicallyActiveAAIDs.append(AAID)
                electronicallyActiveCGSites = copy.deepcopy(electronicallyActiveAAIDs)
                # Now work out what the species is:
                for species, rigidBodies in parameterDict['AARigidBodySpecies'].items():
                    if AAMorphologyDict['body'][electronicallyActiveCGSites[0]] in rigidBodies:
                        self.species = species
                        break
            else:
                raise SystemError("Multiple electronic species defined, but no way to map them without a coarse-grained morphology (no CG morph has been given)")
            self.AAIDs = chromophoreCGSites
        # The position of the chromophore can be calculated easily. Note that here, the `self.image' is the periodic image that the unwrapped_position of the chromophore is located in, relative to the original simulation volume.
        electronicallyActiveUnwrappedPosns = [AAMorphologyDict['unwrapped_position'][AAID] for AAID in electronicallyActiveAAIDs]
        electronicallyActiveTypes = [AAMorphologyDict['type'][AAID] for AAID in electronicallyActiveAAIDs]
        self.unwrappedPosn, self.posn, self.image = self.obtainChromophoreCoM(electronicallyActiveUnwrappedPosns, electronicallyActiveTypes, simDims)
        # A list of the important bonds for this chromophore from the morphology would be useful when determining
        # if a terminating group is already present on this monomer
        self.bonds = self.getImportantBonds(AAMorphologyDict['bond'])
        if CGMorphologyDict is not None:
            # Determine if this chromophore is a repeat unit and therefore will need terminating before ORCA
            CGTypes = set([CGMorphologyDict['type'][CGID] for CGID in chromophoreCGSites])
            # self.terminate = True if any of the CGTypes in this chromophore are defined as having termination conditions in the parameter file
            self.terminate = any(CGType in CGTypes for CGType in [connection[0] for connection in parameterDict['moleculeTerminatingConnections']])
        else:
            try:
                if len(parameterDict['moleculeTerminatingConnections'].keys()) == 0:
                    # Small molecules in atomistic morphology therefore no terminations needed
                    self.terminate = False
            except AttributeError:
                if len(parameterDict['moleculeTerminatingConnections']) == 0:
                    self.terminate = False
            else:
                # No CG morphology, but terminations have been specified, so we're dealing with a polymer
                AATypes = set([AAMorphologyDict['type'][AAID] for AAID in self.AAIDs])
                self.terminate = any(AAType in AATypes for AAType in [connection for connection in parameterDict['moleculeTerminatingConnections']])
        # Now to create a load of placeholder parameters to update later when we have the full list/energy levels
        # The self.neighbours list contains one element for each chromophore within parameterDict['maximumHopDistance']
        # of this one (including periodic boundary conditions). Its format is [[neighbour1ID, relativeImageOfNeighbour1],...]
        self.neighbours = []
        self.dissociationNeighbours = []
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
        pass

    def getImportantBonds(self, bondList):
        importantBonds = []
        for bond in bondList:
            if (bond[1] in self.AAIDs) and (bond[2] in self.AAIDs):
                importantBonds.append(bond)
        return importantBonds

    def obtainChromophoreCoM(self, electronicallyActiveUnwrappedPosns, electronicallyActiveTypes, simDims):
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
    print("Determining chromophores in the system...")
    bondedAtoms = helperFunctions.obtainBondedList(CGMorphologyDict['bond'])
    chromophoreList = [i for i in range(len(CGMorphologyDict['type']))]
    for CGSiteID, chromophoreID in enumerate(chromophoreList):
        CGSiteType = CGMorphologyDict['type'][CGSiteID]
        typesInThisChromophore = [CGSiteType]
        chromophoreList, typesInThisChromophore = updateChromophores(CGSiteID, chromophoreList, bondedAtoms, CGMorphologyDict['type'], typesInThisChromophore, parameterDict)
    chromophoreData = {}
    for atomID, chromoID in enumerate(chromophoreList):
        if chromoID not in list(chromophoreData.keys()):
            chromophoreData[chromoID] = [atomID]
        else:
            chromophoreData[chromoID].append(atomID)
    # Now rename the chromophore IDs so that they increment sensibly (they will be used later
    # for the ORCA files)
    oldKeys = sorted(chromophoreData.keys())
    for newKey, oldKey in enumerate(oldKeys):
        chromophoreData[newKey] = chromophoreData.pop(oldKey)
    print(str(len(list(chromophoreData.keys()))) + " chromophores successfully identified!")
    # Now let's create a list of all the chromophore instances which contain all of the
    # information we could ever want about them.
    chromophoreInstances = []
    for chromoID, chromophoreCGSites in chromophoreData.items():
        print("\rCalculating properties of chromophore %04d of %04d..." % (chromoID, len(list(chromophoreData.keys())) - 1), end=' ')
        sys.stdout.flush()
        chromophoreInstances.append(chromophore(chromoID, chromophoreCGSites, CGMorphologyDict, AAMorphologyDict, CGToAAIDMaster, parameterDict, simDims))
    print("")
    return chromophoreInstances


def calculateChromophoresAA(CGMorphologyDict, AAMorphologyDict, CGToAAIDMaster, parameterDict, simDims, rigidBodies = None):
    # If rigidBodies == None:
    # This function works in the same way as the coarse-grained version above, except 
    # this one iterates through the AA bonds instead. This is FAR SLOWER and so shouldn't
    # be done, except in the case where the coarse-grained morphology does not exist 
    # (because we started with an atomistic morphology and are only interested in running
    # KMC on it)
    # If rigidBodies == AAMorphologyDict['body']:
    # This function uses the rigid bodies specified in parameterDict['AARigidBodySpecies'],
    # and those which have not been specified by iterating through the AA bond list, to determine
    # the chromophores in the system. This is the slowest way to calculate chromophores, but is
    # useful for systems such as BDT-TPD, where there are multiple chromophores of differing
    # species present in the same molecule. As above, this code will only run if an atomistic
    # morphology has been input to MorphCT. If it is coarse-grained, the CG-based "calculateChromophore"
    # function will be used, and will also be a lot faster.
    # The parameterDict['AARigidBodySpecies'] is a dictionary with two keys, 'donor' or 'acceptor'.
    # Each element in the value list corresponds to a new chromophore. These aren't the only atoms
    # that belong to this chromophore, however - there might be a bunch of aliphatic/flexible
    # atoms that are connected, so we need to make sure that we add those too.
    print("Determining chromophores in the system...")
    bondedAtoms = helperFunctions.obtainBondedList(AAMorphologyDict['bond'])
    chromophoreList = [i for i in range(len(AAMorphologyDict['type']))]
    for AASiteID, chromophoreID in enumerate(chromophoreList):
        AASiteType = AAMorphologyDict['type'][AASiteID]
        chromophoreList = updateChromophoresAA(AASiteID, chromophoreList, bondedAtoms, parameterDict, rigidBodies)
    chromophoreData = {}
    for atomID, chromoID in enumerate(chromophoreList):
        if chromoID not in list(chromophoreData.keys()):
            chromophoreData[chromoID] = [atomID]
        else:
            chromophoreData[chromoID].append(atomID)
    # Now rename the chromophore IDs so that they increment sensibly (they will be used later
    # for the ORCA files)
    oldKeys = sorted(chromophoreData.keys())
    for newKey, oldKey in enumerate(oldKeys):
        chromophoreData[newKey] = chromophoreData.pop(oldKey)
    print(str(len(list(chromophoreData.keys()))) + " chromophores successfully identified!")
    # Now let's create a list of all the chromophore instances which contain all of the
    # information we could ever want about them.
    chromophoreInstances = []
    for chromoID, chromophoreCGSites in chromophoreData.items():
        print("\rCalculating properties of chromophore %04d of %04d..." % (chromoID, len(list(chromophoreData.keys())) - 1), end=' ')
        sys.stdout.flush()
        chromophoreInstances.append(chromophore(chromoID, chromophoreCGSites, CGMorphologyDict, AAMorphologyDict, CGToAAIDMaster, parameterDict, simDims))
    print("")
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


def updateChromophoresAA(atomID, chromophoreList, bondedAtoms, parameterDict, rigidBodies = None):
    # This version of the update chromophores function does not check for CG site types, instead
    # just adding all bonded atoms. Therefore it should only be used in the case of already-atomistic
    # morphologies (no CG morph specified) containing ONLY small molecules
    try:
        for bondedAtom in bondedAtoms[atomID]:
            if rigidBodies is not None:
                # Skip if the bonded atom belongs to a different rigid body
                if ((rigidBodies[bondedAtom] != -1) and (rigidBodies[atomID] != -1)) and (rigidBodies[bondedAtom] != rigidBodies[atomID]):
                    continue
            # If the atomID of the bonded atom is larger than that of the current one,
            # update the bonded atom's ID to the current one's to put it in this chromophore,
            # then iterate through all of the bonded atom's neighbours
            if chromophoreList[bondedAtom] > chromophoreList[atomID]:
                chromophoreList[bondedAtom] = chromophoreList[atomID]
                chromophoreList = updateChromophoresAA(bondedAtom, chromophoreList, bondedAtoms, parameterDict, rigidBodies)
            # If the atomID of the current atom is larger than that of the bonded one,
            # update the current atom's ID to the bonded one's to put it in this chromophore,
            # then iterate through all of the current atom's neighbours
            elif chromophoreList[bondedAtom] < chromophoreList[atomID]:
                chromophoreList[atomID] = chromophoreList[bondedAtom]
                chromophoreList = updateChromophoresAA(atomID, chromophoreList, bondedAtoms, parameterDict, rigidBodies)
            # Else: both the current and the bonded atom are already known to be in this
            # chromophore, so we don't have to do anything else.
    except KeyError:
        # This means that there are no bonded CG sites (i.e. it's a single chromophore)
        pass
    return chromophoreList


def determineNeighbours(chromophoreList, parameterDict, simDims):
    for chromophore1 in chromophoreList:
        print("\rIdentifying neighbours of chromophore %04d of %04d..." % (chromophore1.ID, len(chromophoreList) - 1), end=' ')
        sys.stdout.flush()
        for chromophore2 in chromophoreList:
            # Skip if chromo2 is chromo1
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
                chromo1DissociationNeighbourIDs = [neighbourData[0] for neighbourData in chromophore1.dissociationNeighbours]
                chromo2DissociationNeighbourIDs = [neighbourData[0] for neighbourData in chromophore2.dissociationNeighbours]
                # Also, make the deltaE and the Tij lists as long as the neighbour lists for easy access later
                if (chromophore1.species == chromophore2.species) and (chromophore2.ID not in chromo1NeighbourIDs):
                    chromophore1.neighbours.append([chromophore2.ID, relativeImageOfChromo2])
                    chromophore1.neighboursDeltaE.append(None)
                    chromophore1.neighboursTI.append(None)
                elif (chromophore1.species != chromophore2.species) and (chromophore2.ID not in chromo2DissociationNeighbourIDs):
                    chromophore1.dissociationNeighbours.append([chromophore2.ID, relativeImageOfChromo2])
                if (chromophore1.species == chromophore2.species) and (chromophore1.ID not in chromo2NeighbourIDs):
                    chromophore2.neighbours.append([chromophore1.ID, list(-np.array(relativeImageOfChromo2))])
                    chromophore2.neighboursDeltaE.append(None)
                    chromophore2.neighboursTI.append(None)
                elif (chromophore1.species != chromophore2.species) and (chromophore1.ID not in chromo1DissociationNeighbourIDs):
                    chromophore2.dissociationNeighbours.append([chromophore1.ID, list(-np.array(relativeImageOfChromo2))])
        # DEBUG TESTING
        # if chromophore1.ID == 1961:
        #     print ""
        #     print chromophore1.posn
        #     print chromophore1.AAIDs
        #     for neighbour in chromophore1.neighbours:
        #         print neighbour[0], neighbour[1], chromophoreList[neighbour[0]].posn
        #     plotChromoNeighbours(chromophore1, chromophoreList, simDims)
    print("")
    return chromophoreList


def execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList):
    simDims = [[-AAMorphologyDict['lx'] / 2.0, AAMorphologyDict['lx'] / 2.0], [-AAMorphologyDict['ly'] / 2.0, AAMorphologyDict['ly'] / 2.0], [-AAMorphologyDict['lz'] / 2.0, AAMorphologyDict['lz'] / 2.0]]
    if len(parameterDict['CGToTemplateDirs']) > 0:
        # Normal operation using the coarse-grained morphology
        chromophoreList = calculateChromophores(CGMorphologyDict, AAMorphologyDict, CGToAAIDMaster, parameterDict, simDims)
    elif (len(parameterDict['CGSiteSpecies']) == 1) and (len(parameterDict['AARigidBodySpecies']) == 0):
        # Small molecule system with only one electronic species
        chromophoreList = calculateChromophoresAA(CGMorphologyDict, AAMorphologyDict, CGToAAIDMaster, parameterDict, simDims)
    else:
        # Other system, with electronically active species specified as rigid bodies using AARigidBodySpecies in parameter file
        chromophoreList = calculateChromophoresAA(CGMorphologyDict, AAMorphologyDict, CGToAAIDMaster, parameterDict, simDims, rigidBodies=AAMorphologyDict['body'])
    chromophoreList = determineNeighbours(chromophoreList, parameterDict, simDims)
    # Now we have updated the chromophoreList, rewrite the pickle with this new information.
    #### SANITY CHECK  ####
    for index, chromo in enumerate(chromophoreList):
        if index != chromo.ID:
            print("Inconsistency found in the ordering of the chromophoreList, rewriting the chromophoreList in the correct order...")
            newChromophoreList = []
            for chromo in chromophoreList:
                newChromophoreList.append(0)
            for chromo in chromophoreList:
                newChromophoreList[chromo.ID] = chromo
            chromophoreList = newChromophoreList
            break
    #### END OF SANITY CHECK ####
    pickleName = parameterDict['outputMorphDir'] + '/' + parameterDict['morphology'][:-4] + '/code/' + parameterDict['morphology'][:-4] + '.pickle'
    helperFunctions.writePickle((AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList), pickleName)
    return AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList


if __name__ == "__main__":
    try:
        pickleFile = sys.argv[1]
    except:
        print("Please specify the pickle file to load to continue the pipeline from this point.")
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = helperFunctions.loadPickle(pickleFile)
    execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList)
