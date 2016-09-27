import copy
import numpy as np
import helperFunctions
import cPickle as pickle


class morphology:
    def __init__(self, morphologyXML, morphologyName, parameterDict):
        self.parameterDict = parameterDict
        for key, value in parameterDict.iteritems():
            self.__dict__[key] = value
        self.xmlPath = morphologyXML
        self.morphologyName = morphologyName
        # self.inputSigma is the `compression value' in Angstroms that has been used to scale the morphology
        # E.G. the P3HT template uses sigma = 1, but the Marsh morphologies use sigma = 3.
        self.CGDictionary = helperFunctions.loadMorphologyXML(self.xmlPath, sigma=self.inputSigma)
        self.CGDictionary = helperFunctions.addUnwrappedPositions(self.CGDictionary)

    def analyseMorphology(self):
        print "Finding molecules..."
        moleculeIDs, moleculeLengths = self.splitMolecules()
        rollingAAIndex = 0
        boxSize = [self.CGDictionary['lx'], self.CGDictionary['ly'], self.CGDictionary['lz']]
        CGMorphologyDict = {}
        AAMorphologyDict = {}
        for boxDimension in ['lx', 'ly', 'lz']:
            CGMorphologyDict[boxDimension] = self.CGDictionary[boxDimension]
            AAMorphologyDict[boxDimension] = self.CGDictionary[boxDimension]
        CGtoAAIDMaster = []  # This is a list of dictionaries. Elements in the list correspond to molecules
        # (useful for splitting out individual molecules for the xyz conversion) within the element, the
        # dictionary key is the CG site, the value is a list containing the CG type (e.g. 'thio') as the
        # first element and then another list of all the AAIDs corresponding to that CG site as the second
        # element.

        # Create a ghost particle dictionary to be added at the end of the morphology.
        # This way, we don't mess up the order of atoms later on when trying to split
        # back up into individual molecules and monomers.
        # The ghost dictionary contains all of the type T and type X particles that will
        # anchor the thiophene rings to the CG COM positions.
        ghostDictionary = {'position': [], 'image': [], 'unwrapped_position': [], 'mass': [], 'diameter': [], 'type': [], 'body': [], 'bond': [], 'angle': [], 'dihedral': [], 'improper': [], 'charge': []}
        print "Adding molecules to the system..."
        for moleculeNumber in range(len(moleculeIDs)):
            print "Adding molecule number", moleculeNumber, "\r",
            # print "Rolling AA Index =", rollingAAIndex
            CGMoleculeDict, AAMoleculeDict, CGtoAAIDs, ghostDictionary = atomistic(moleculeNumber, moleculeIDs[moleculeNumber], self.CGDictionary, moleculeLengths, rollingAAIndex, ghostDictionary, self.parameterDict).returnData()
            CGtoAAIDMaster.append(CGtoAAIDs)
            for key in CGMoleculeDict.keys():
                if key not in ['lx', 'ly', 'lz']:
                    if key not in CGMorphologyDict.keys():
                        CGMorphologyDict[key] = CGMoleculeDict[key]
                    else:
                        CGMorphologyDict[key] += CGMoleculeDict[key]
                    if key not in AAMorphologyDict.keys():
                        AAMorphologyDict[key] = AAMoleculeDict[key]
                    else:
                        AAMorphologyDict[key] += AAMoleculeDict[key]
            rollingAAIndex += len(AAMoleculeDict['type'])
        # Now add the ghost dictionary to the end of the morphology file
        totalNumberOfAtoms = len(AAMorphologyDict['type'])  # Should be == rollingAAIndex, but don't want to take any chances
        # Add in the wrapped positions of the ghosts. Need to know sim dims for this
        for key in ['lx', 'ly', 'lz']:
            ghostDictionary[key] = AAMorphologyDict[key]
        ghostDictionary = helperFunctions.addWrappedPositions(ghostDictionary)
        for key in ['lx', 'ly', 'lz']:
            ghostDictionary.pop(key)
        # The real atoms that the ghost particles are bonded to are already correct and no longer need to be changed.
        # However, the ghost particles themselves will have the wrong indices if we were to add them to the system directly.
        # Therefore, increment all of the ghost bond indices that begin with a * (ghost particle) by the total number of atoms already in the system.
        for bondNo, bond in enumerate(ghostDictionary['bond']):
            if str(bond[1])[0] == '*':
                ghostDictionary['bond'][bondNo][1] = int(bond[1][1:]) + totalNumberOfAtoms
            if str(bond[2])[0] == '*':
                ghostDictionary['bond'][bondNo][2] = int(bond[2][1:]) + totalNumberOfAtoms
        # Now append all ghosts to morphology
        for key in ghostDictionary.keys():
            AAMorphologyDict[key] += ghostDictionary[key]
        # Finally, update the number of atoms
        AAMorphologyDict['natoms'] += len(ghostDictionary['type'])
        print "\n"
        print "Writing XML file..."
        AAFileName = './outputFiles/' + self.morphologyName + '/morphology/' + self.morphologyName + '.xml'
        writeXML(AAMorphologyDict, './templates/template.xml', AAFileName)
        toPickle = (AAFileName, CGMorphologyDict, AAMorphologyDict, CGtoAAIDMaster, [], boxSize)
        print "Writing pickle file..."
        pickleFileName = './outputFiles/' + self.morphologyName + '/morphology/' + self.morphologyName + '.pickle'
        with open(pickleFileName, 'w+') as pickleFile:
            pickle.dump(toPickle, pickleFile)
        print "Pickle file written to", pickleFileName
        return AAFileName, CGMorphologyDict, AAMorphologyDict, CGtoAAIDMaster, [], boxSize
        # ### Before we split into segments, let's fine grain the molecules (and write an xml so that we can see the output in VMD)
        # rollingMoleculeNumber = 0
        # for moleculeAtoms in moleculeIDs:
        #     moleculeNumber = str(rollingMoleculeNumber)
        #     while len(moleculeNumber) < 4:
        #         moleculeNumber = '0'+moleculeNumber
        #     # Check to see if pickle file is already present, skip if it is
        #     pickleFound = False
        #     try:
        #         currentMolDirContents = os.listdir('./outputChains/'+morphologyName+'/mol'+moleculeNumber)
        #         for fileName in currentMolDirContents:
        #             if '.pickle' in fileName:
        #                 # Pickle present so skip
        #                 pickleFound = True
        #                 break
        #         if pickleFound == True:
        #             rollingMoleculeNumber += 1
        #             # GET ONE MOLECULE WORKING COMPLETELY FIRST
        #             # CONTINUE TO NEXT MOLECULE
        #             continue
        #     except OSError:
        #         # Directory doesn't exist
        #         pass
        #     AAfileName, CGMoleculeDict, AAMoleculeDict, CGtoAAIDs, boxSize = atomistic(moleculeAtoms, self.CGDictionary, morphologyName, moleculeNumber).returnData()
        #     toPickle = (AAfileName, CGMoleculeDict, CGtoAAIDs, boxSize)
        #     pickleFileName = './outputChains/'+morphologyName+'/mol'+moleculeNumber+'/mol'+moleculeNumber+'.pickle'
        #     with open(pickleFileName, 'w+') as pickleFile:
        #         pickle.dump(toPickle, pickleFile)
        #     # print "~~~ THIS WILL EVENTUALLY LOOP OVER ALL MOLECULES IN THE MORPHOLOGY, GENERATING THE REQUIRED XML AND PICKLE FILES ~~~"
        #     # print "~~~ FOR NOW JUST GET ONE MOLECULE WORKING IN THE PIPELINE FIRST ~~~"
        #     # runHoomd.execute(AAfileName, CGMoleculeDict, CGtoAAIDs, boxSize)
        #     rollingMoleculeNumber += 1
        #     # JUST RUN THE FIRST 10 MOLECULES TO CHECK EVERYTHING IS OK
        #     # if rollingMoleculeNumber == 9:
        #     #     return 0
        #     # raise SystemError('GET ONE MOLECULE WORKING COMPLETELY FIRST')

        #####
        # morphologySegments = []
        # for moleculeAtoms in moleculeIDs:
        #    print "WHEN CALLING MOLECULE CLASS, INPUT THE MOLECULE ENDS FROM THE CALCULATION IN THE ATOMISTIC CLASS TO SAVE ON CODE LINES"
        #    morphologySegments.append(molecule(moleculeAtoms, self.CGDictionary).returnSegments())
        # CALCULATE AVERAGE SEGMENT LENGTH FOR THE MORPHOLOGY
        # totalSegments = 0
        # segmentLength = 0
        # for moleculeSegments in morphologySegments:
        #     for segment in moleculeSegments:
        #         totalSegments += 1
        #         segmentLength += len(segment)
        # print "Average Segment Length =", segmentLength/float(totalSegments)

    def splitMolecules(self):
        moleculeList = []
        moleculeLengths = []
        bondList = copy.deepcopy(self.CGDictionary['bond'])
        while len(bondList) > 0:
            # Add the first two CG sites of the molecule
            thisMolecule = [bondList[0][1], bondList[0][2]]
            addedNewSite = True
            currentNumberOfSitesInMolecule = 2
            while addedNewSite is True:
                addedNewSite = False
                bondPopList = []
                for bondNo, bond in enumerate(bondList):
                    if (bond[1] in thisMolecule) and (bond[2] not in thisMolecule):
                        thisMolecule.append(bond[2])
                        currentNumberOfSitesInMolecule += 1
                        addedNewSite = True
                    elif (bond[2] in thisMolecule) and (bond[1] not in thisMolecule):
                        thisMolecule.append(bond[1])
                        currentNumberOfSitesInMolecule += 1
                        addedNewSite = True
                    elif (bond[1] in thisMolecule) and (bond[2] in thisMolecule):
                        pass
                    else:
                        continue
                    bondPopList.append(bondNo)
                bondPopList.sort(reverse=True)
                for bondIndex in bondPopList:
                    bondList.pop(bondIndex)
            moleculeList.append(thisMolecule)
            # Length of the molecule is the total number of CG sites / number of CG sites in each repeat unit
            # Sanity check:
            if currentNumberOfSitesInMolecule % len(self.CGToTemplateAAIDs.keys()) != 0:
                raise SystemError("Issue with splitting the morphology into molecules - does the template file contain exactly one repeat unit (minus any terminating groups)?")
            moleculeLength = np.round(currentNumberOfSitesInMolecule / len(self.CGToTemplateAAIDs.keys()))  # Just in case of floating point errors
            moleculeLengths.append(moleculeLength)
        return moleculeList, moleculeLengths


class molecule:
    def __init__(self, moleculeIDs, CGDictionary):
        # TOLERANCE SHOULD BE IN A PARAMETER FILE #
        self.tolerance = np.pi / 6.
        #################################################
        self.atomIDs = moleculeIDs
        self.CGDictionary = CGDictionary
        print "Identifying Segments..."
        self.segments = self.findSegments()

    def returnSegments(self):
        return self.segments

    def findSegments(self):
        segmentMaster = [[]]
        polymerBackboneIDs, moleculeEnds = self.obtainBackboneAtoms()
        atomUnderConsideration = moleculeEnds[0]
        segmentMaster[-1].append(atomUnderConsideration)
        firstAtomInSegment = True
        while True:
            print "\nCurrent atom under consideration =", atomUnderConsideration, "Molecule ends =", moleculeEnds, "New Segment =", firstAtomInSegment
            if atomUnderConsideration == moleculeEnds[1]:
                # Continue until the end of the current molecule
                break
            atomUnderConsiderationCoords = self.CGDictionary['unwrapped_position'][atomUnderConsideration]
            neighbouringBackboneAtoms = self.findBondedNeighbours(atomUnderConsideration, polymerBackboneIDs)
            for atom in neighbouringBackboneAtoms:
                # Should give next atom along in the chain
                if atom in segmentMaster[-1]:
                    # Atom already in current segment
                    continue
                neighbouringBackboneAtom = atom
            neighbouringBackboneAtomCoords = self.CGDictionary['unwrapped_position'][neighbouringBackboneAtom]
            print "Current Atom Posn =", self.CGDictionary['position'][atomUnderConsideration], self.CGDictionary['image'][atomUnderConsideration]
            print "Neighbour Atom Posn =", self.CGDictionary['position'][atomUnderConsideration], self.CGDictionary['image'][atomUnderConsideration]
            if firstAtomInSegment is True:
                axisVector = helperFunctions.findAxis(atomUnderConsiderationCoords, neighbouringBackboneAtomCoords)
                firstAtomInSegment = False
                print "No longer first atom in segment"
            separationVector = helperFunctions.findAxis(atomUnderConsiderationCoords, neighbouringBackboneAtomCoords)
            dotProduct = np.dot(axisVector, separationVector)
            print dotProduct
            # print "Axis Vector for this atom =", axisVector, "Separation Vector for this atom =", separationVector, "Dot Product for this atom =", dotProduct
            if abs(dotProduct - 1.0) <= 1E-8:
                # Floating point inaccuracy check
                dotProduct = 1.0
            separationAngle = np.arccos(dotProduct)
            print "Pos 1", atomUnderConsiderationCoords, "Pos 2", neighbouringBackboneAtomCoords
            if abs(separationAngle) <= self.tolerance:
                # Atom belongs in this segment and coherence hasn't been disrupted
                pass
            else:
                # Orbital conjugation has been disrupted, end segment and create a new one
                print "Separation Angle =", separationAngle, ">", self.tolerance, "therefore time for a new segment."
                segmentMaster.append([])
                firstAtomInSegment = True
            segmentMaster[-1].append(neighbouringBackboneAtom)
            atomUnderConsideration = neighbouringBackboneAtom
            # The new axis Vector becomes the previous backbone vector to be used in the next loop iteration
            axisVector = helperFunctions.findAxis(atomUnderConsiderationCoords, neighbouringBackboneAtomCoords)
        return segmentMaster

    def findBondedNeighbours(self, atomUnderConsideration, backboneAtoms):
        bondedNeighbours = []
        for bond in self.CGDictionary['bond']:
            if bond[1] == atomUnderConsideration:
                if bond[2] in backboneAtoms:
                    bondedNeighbours.append(bond[2])
            elif bond[2] == atomUnderConsideration:
                if bond[1] in backboneAtoms:
                    bondedNeighbours.append(bond[1])
        return bondedNeighbours

    def obtainBackboneAtoms(self):
        polymerBackboneIDs = []
        for atomID in self.atomIDs:
            if self.CGDictionary['type'][atomID] == 'A':
                polymerBackboneIDs.append(atomID)
        # Find the ends:
        moleculeEnds = []
        numberOfBonds = {}
        for atom in polymerBackboneIDs:
            numberOfBonds[atom] = 0
        for bond in self.CGDictionary['bond']:
            if bond[0] == 'bondA':
                if (bond[1] in numberOfBonds) and (bond[2] in numberOfBonds):
                    numberOfBonds[bond[1]] += 1
                    numberOfBonds[bond[2]] += 1
        for atomID, bondQuantity in numberOfBonds.iteritems():
            if bondQuantity == 1:
                moleculeEnds.append(atomID)
        moleculeEnds.sort()
        return polymerBackboneIDs, moleculeEnds


class atomistic:
    def __init__(self, moleculeIndex, siteIDs, CGDictionary, moleculeLengths, rollingAAIndex, ghostDictionary, parameterDict):
        # This class sees individual molecules.
        self.noAtomsInMorphology = rollingAAIndex
        self.moleculeIndex = moleculeIndex
        self.moleculeLengths = moleculeLengths
        self.siteIDs = siteIDs
        self.CGDictionary = CGDictionary
        self.CGMonomerDictionary = self.getCGMonomerDict()
        for key, value in parameterDict.iteritems():
            self.__dict__[key] = value
        AATemplateDictionary = helperFunctions.loadMorphologyXML(self.AATemplateFile)
        self.AATemplateDictionary = helperFunctions.addUnwrappedPositions(AATemplateDictionary)
        # thioAA, alk1AA, alk2AA are now all unused. Instead, use the CGToTemplateAAIDs dictionary which can have arbitrary length
        self.AADictionary, self.atomIDLookupTable, self.ghostDictionary = self.runFineGrainer(ghostDictionary)

    def returnData(self):
        return self.CGMonomerDictionary, self.AADictionary, self.atomIDLookupTable, self.ghostDictionary

    def obtainBackboneAtoms(self):
        polymerBackboneIDs = []
        for atomID in self.atomIDs:
            if self.CGDictionary['type'][atomID] == 'A':
                polymerBackboneIDs.append(atomID)
        # Find the ends:
        moleculeEnds = []
        numberOfBonds = {}
        for atom in polymerBackboneIDs:
            numberOfBonds[atom] = 0
        for bond in self.CGDictionary['bond']:
            if bond[0] == 'bondA':
                if (bond[1] in numberOfBonds) and (bond[2] in numberOfBonds):
                    numberOfBonds[bond[1]] += 1
                    numberOfBonds[bond[2]] += 1
        for atomID, bondQuantity in numberOfBonds.iteritems():
            if bondQuantity == 1:
                moleculeEnds.append(atomID)
        moleculeEnds.sort()
        return polymerBackboneIDs, moleculeEnds

    def getCGMonomerDict(self):
        CGMonomerDictionary = {'position': [], 'image': [], 'mass': [], 'diameter': [], 'type': [], 'body': [], 'bond': [], 'angle': [], 'dihedral': [], 'improper': [], 'charge': [], 'lx': 0, 'ly': 0, 'lz': 0}
        # First, do just the positions and find the newsiteIDs for each CG site
        for siteID in self.siteIDs:
            CGMonomerDictionary['position'].append(self.CGDictionary['position'][siteID])
        # Now sort out the other one-per-atom properties
        for key in ['image', 'mass', 'diameter', 'type', 'body', 'charge']:
            if len(self.CGDictionary[key]) != 0:
                for siteID in self.siteIDs:
                    CGMonomerDictionary[key].append(self.CGDictionary[key][siteID])
        # Now rewrite the bonds based on the newsiteIDs
        for key in ['bond', 'angle', 'dihedral', 'improper']:
            for element in self.CGDictionary[key]:
                for siteID in self.siteIDs:
                    if (siteID in element) and (element not in CGMonomerDictionary[key]):
                        CGMonomerDictionary[key].append(element)
        # Now update the box parameters
        for key in ['lx', 'ly', 'lz']:
            CGMonomerDictionary[key] = self.CGDictionary[key]
        CGMonomerDictionary = helperFunctions.addUnwrappedPositions(CGMonomerDictionary)
        CGMonomerDictionary['natoms'] = len(CGMonomerDictionary['position'])
        return CGMonomerDictionary

    def runFineGrainer(self, ghostDictionary):
        AADictionary = {'position': [], 'image': [], 'unwrapped_position': [], 'mass': [], 'diameter': [], 'type': [], 'body': [], 'bond': [], 'angle': [], 'dihedral': [], 'improper': [], 'charge': [], 'lx': 0, 'ly': 0, 'lz': 0}
        CGCoMs = self.getAATemplatePosition(self.CGToTemplateAAIDs)
        # Need to keep track of the atom ID numbers globally - runFineGrainer sees individual monomers, atomistic sees molecules and the XML needs to contain the entire morphology.
        noAtomsInMolecule = 0
        CGTypeList = {}
        for siteID in self.siteIDs:
            CGTypeList[siteID] = self.CGDictionary['type'][siteID]
        monomerList = self.sortIntoMonomers(CGTypeList)
        currentMonomerIndex = sum(self.moleculeLengths[:self.moleculeIndex])
        atomIDLookupTable = {}
        totalPermittedAtoms = len(self.AATemplateDictionary['type']) * len(monomerList)
        for monomerNo, monomer in enumerate(monomerList):
            thisMonomerDictionary = copy.deepcopy(self.AATemplateDictionary)
            for key in ['lx', 'ly', 'lz']:
                thisMonomerDictionary[key] = self.CGDictionary[key]
            if len(thisMonomerDictionary['image']) == 0:
                thisMonomerDictionary['image'] = [[0, 0, 0]] * len(thisMonomerDictionary['position'])
            for siteID in monomer:
                sitePosn = np.array(self.CGDictionary['unwrapped_position'][siteID])
                siteTranslation = sitePosn - CGCoMs[CGTypeList[siteID]]
                atomIDLookupTable[siteID] = [CGTypeList[siteID], [x + noAtomsInMolecule + self.noAtomsInMorphology for x in self.CGToTemplateAAIDs[CGTypeList[siteID]]]]
                # Add the atoms in based on the CG site position
                for AAID in self.CGToTemplateAAIDs[CGTypeList[siteID]]:
                    thisMonomerDictionary['unwrapped_position'][AAID] = list(np.array(thisMonomerDictionary['unwrapped_position'][AAID]) + siteTranslation)
                # Next sort out the rigid bodies
                if CGTypeList[siteID] in self.rigidBodySites:
                    # Every rigid body needs a ghost particle that describes its CoM
                    AAIDPositions = []
                    AAIDAtomTypes = []
                    for AAID in self.rigidBodySites[CGTypeList[siteID]]:
                        thisMonomerDictionary['body'][AAID] = currentMonomerIndex
                        AAIDPositions.append(thisMonomerDictionary['unwrapped_position'][AAID])
                        AAIDAtomTypes.append(thisMonomerDictionary['type'][AAID])
                    # Now create the ghost particle describing the rigid body
                    ghostDictionary['unwrapped_position'].append(helperFunctions.calcCOM(AAIDPositions, listOfAtomTypes=AAIDAtomTypes))
                    ghostDictionary['mass'].append(1.0)
                    ghostDictionary['diameter'].append(1.0)
                    ghostDictionary['type'].append('R' + str(CGTypeList[siteID]))
                    ghostDictionary['body'].append(currentMonomerIndex)
                    ghostDictionary['charge'].append(0.0)
                    # Then create the corresponding CG anchorpoint
                    ghostDictionary['unwrapped_position'].append(self.CGDictionary['unwrapped_position'][siteID])
                    ghostDictionary['mass'].append(1.0)
                    ghostDictionary['diameter'].append(1.0)
                    ghostDictionary['type'].append('X' + str(CGTypeList[siteID]))
                    ghostDictionary['body'].append(-1)
                    ghostDictionary['charge'].append(0.0)
                    # Now create a bond between them
                    # We want to bond together the previous two ghost particles, so this should work as it requires no knowledge of the number of ghost particles already in the system.
                    ghostDictionary['bond'].append([str(ghostDictionary['type'][-2]) + '-' + str(ghostDictionary['type'][-1]), '*' + str(len(ghostDictionary['type']) - 2), '*' + str(len(ghostDictionary['type']) - 1)])
                else:
                    # Create a ghost particle that describe the CG anchorpoint for the non-rigid body
                    ghostDictionary['unwrapped_position'].append(self.CGDictionary['unwrapped_position'][siteID])
                    ghostDictionary['mass'].append(1.0)
                    ghostDictionary['diameter'].append(1.0)
                    ghostDictionary['type'].append('X' + str(CGTypeList[siteID]))
                    ghostDictionary['body'].append(-1)
                    ghostDictionary['charge'].append(0.0)
                    # Add in bonds between the CG anchorpoints and the atom closest to the CG site
                    # Find the atom closest to the CG site
                    closestAtomID = None
                    closestAtomPosn = 1E99
                    for AAID, AAPosition in enumerate(thisMonomerDictionary['unwrapped_position']):
                        separation = helperFunctions.calculateSeparation(self.CGDictionary['unwrapped_position'][siteID], AAPosition)
                        if separation < closestAtomPosn:
                            closestAtomPosn = separation
                            closestAtomID = AAID
                    # Add in the bond:
                    # Note that, in order to distinguish between the ghostAtomIDs and the realAtomIDs, I've put an underscore in front of the closestAtomID, and a * in front of the ghostAtomID.
                    # When incrementing the atomIDs for this monomer or this molecule, the readlAtomIDs will be incremented correctly.
                    # Later, when the ghost atoms are added to the main system, the ghostAtomIDs will be incremented according to the number of atoms in the whole system (i.e. the ghosts appear at the end of the real atoms).
                    # At this time, the realAtomIDs will be left unchanged because they are already correct for the whole system.
                    ghostDictionary['bond'].append([str(ghostDictionary['type'][-1]) + '-' + str(thisMonomerDictionary['type'][closestAtomID]), '*' + str(len(ghostDictionary['type']) - 1), '_' + str(closestAtomID + noAtomsInMolecule)])
            thisMonomerDictionary = helperFunctions.addWrappedPositions(thisMonomerDictionary)
            thisMonomerDictionary = helperFunctions.addMasses(thisMonomerDictionary)
            thisMonomerDictionary = helperFunctions.addDiameters(thisMonomerDictionary)
            # Now add in the bonds between CGSites in this monomer
            for bond in self.CGDictionary['bond']:
                if (bond[1] in monomer) and (bond[2] in monomer):
                    CGBondType = bond[0]
                    thisMonomerDictionary['bond'].append(self.CGToTemplateBonds[CGBondType])
            # Now need to add in the additionalConstraints for this monomer (which include the bond, angle and dihedral for the inter-monomer connections.
            # However, we need a check to make sure that we don't add stuff for the final monomer (because those atoms +25 don't exist in this molecule!)
            for constraint in self.additionalConstraints:
                # Check that we're not at the final monomer
                atFinalMonomer = False
                for atomID in constraint[1:]:
                    if (noAtomsInMolecule + atomID + 1) > totalPermittedAtoms:
                        atFinalMonomer = True
                        break
                if atFinalMonomer is True:
                    break
                # Work out which key to write the constraint to based on its length:
                # 3 = Bond, 4 = Angle, 5 = Dihedral, 6 = Improper
                constraintType = ['bond', 'angle', 'dihedral', 'improper']
                thisMonomerDictionary[constraintType[len(constraint) - 3]].append(constraint)
            # Finally, increment the atom IDs to take into account previous monomers in this molecule and then update the AADictionary.
            # Note that the ghost dictionary bond was already updated to have the correct realAtom AAID for this molecule when the bond was created. Therefore, leave the ghost dictionary unchanged.
            thisMonomerDictionary, ghostDictionary = helperFunctions.incrementAtomIDs(thisMonomerDictionary, ghostDictionary, noAtomsInMolecule, modifyGhostDictionary=False)
            # Find the connecting atoms to the terminating units based on monomer number
            if len(self.moleculeTerminatingConnections) != 0:
                if monomerNo == 0:
                    startAtomIndex = noAtomsInMolecule + self.moleculeTerminatingConnections[0][1]
                elif monomerNo == len(monomerList) - 1:
                    endAtomIndex = noAtomsInMolecule + self.moleculeTerminatingConnections[1][1]
            noAtomsInMolecule += len(thisMonomerDictionary['type'])
            currentMonomerIndex += 1
            AADictionary = self.updateMoleculeDictionary(thisMonomerDictionary, AADictionary)
        # All Monomers sorted, now for the final bits
        AADictionary['natoms'] = noAtomsInMolecule
        for key in ['lx', 'ly', 'lz']:
            AADictionary[key] = thisMonomerDictionary[key]
        # Add in the terminating units
        # TODO: This code only permits hydrogens to be used as terminating units current. Perhaps it wouldn't be that hard to implement a template-based termination unit for enhanced flexibility.
        if len(self.moleculeTerminatingConnections) != 0:
            AADictionary, startTerminatingHydrogen = helperFunctions.addTerminatingHydrogen(AADictionary, startAtomIndex)
            AADictionary, endTerminatingHydrogen = helperFunctions.addTerminatingHydrogen(AADictionary, endAtomIndex)
            atomIDLookupTable[monomerList[0][0]].append(startTerminatingHydrogen + self.noAtomsInMorphology)
            atomIDLookupTable[monomerList[-1][0]].append(endTerminatingHydrogen + self.noAtomsInMorphology)
        # Now the molecule is done, we need to add on the correct identifying numbers for all the bonds, angles and dihedrals
        # (just as we did between monomers) for the other molecules in the system, so that they all connect to the right atoms
        # Note that here we need to increment the '_'+ATOMIDs in the ghost dictionary to take into account the number of molecules.
        AADictionary, ghostDictionary = helperFunctions.incrementAtomIDs(AADictionary, ghostDictionary, self.noAtomsInMorphology, modifyGhostDictionary=True)
        return AADictionary, atomIDLookupTable, ghostDictionary

    def sortIntoMonomers(self, typeListSequence):
        monomerList = []
        moleculeSiteIDs = copy.deepcopy(self.siteIDs)
        bondList = copy.deepcopy(self.CGDictionary['bond'])
        while len(moleculeSiteIDs) > 0:
            # Add the first atom to the first monomer
            thisMonomer = []
            monomerStartSiteType = typeListSequence[moleculeSiteIDs[0]]
            thisMonomer.append(moleculeSiteIDs[0])
            addedNewSite = True
            while addedNewSite is True:
                addedNewSite = False
                bondPopList = []
                # Find bonded atoms that are not of the same type
                for bondNo, bond in enumerate(bondList):
                    if (bond[1] in thisMonomer) and (bond[2] not in thisMonomer) and (typeListSequence[bond[2]] != monomerStartSiteType):
                        thisMonomer.append(bond[2])
                        addedNewSite = True
                    elif (bond[2] in thisMonomer) and (bond[1] not in thisMonomer) and (typeListSequence[bond[1]] != monomerStartSiteType):
                        thisMonomer.append(bond[1])
                        addedNewSite = True
                    elif (bond[1] in thisMonomer) and (bond[2] in thisMonomer):
                        pass
                    else:
                        continue
                    bondPopList.append(bondNo)
                bondPopList.sort(reverse=True)
                for bondIndex in bondPopList:
                    bondList.pop(bondIndex)
            for atomIndex in thisMonomer:
                moleculeSiteIDs.remove(atomIndex)
            monomerList.append(thisMonomer)
        monomerTypesList = []
        for monomer in monomerList:
            monomerTypesList.append([])
            for atom in monomer:
                monomerTypesList[-1].append(typeListSequence[atom])
        return monomerList

    def updateMoleculeDictionary(self, currentMonomerDictionary, AADictionary):
        keyList = AADictionary.keys()
        keyList.remove('lx')
        keyList.remove('ly')
        keyList.remove('lz')
        for key in keyList:
            for value in currentMonomerDictionary[key]:
                AADictionary[key].append(value)
        return AADictionary

    def getAATemplatePosition(self, CGToTemplateAAIDs):
        CGCoMs = {}
        for siteName in CGToTemplateAAIDs.keys():
            atomIDs = CGToTemplateAAIDs[siteName]
            sitePositions = []
            siteTypes = []
            for atomID in atomIDs:
                siteTypes.append(self.AATemplateDictionary['type'][atomID])
                sitePositions.append(self.AATemplateDictionary['unwrapped_position'][atomID])
            # These output as numpy arrays because we can't do maths with lists
            CGCoMs[siteName] = helperFunctions.calcCOM(sitePositions, listOfAtomTypes=siteTypes)
        return CGCoMs


class writeXML:
    def __init__(self, dataDictionary, templateFile, outputFile):
        self.writeUnwrappedPositionsOnly = True
        self.dataDictionary = dataDictionary
        self.templateFile = templateFile
        self.outputFile = outputFile
        xmlTemplateData = self.loadTemplate()
        newTemplateData = self.updateTemplate(xmlTemplateData)
        self.writeData(newTemplateData)

    def loadTemplate(self):
        with open(self.templateFile, 'r') as xmlFile:
            xmlTemplateData = xmlFile.readlines()
        return xmlTemplateData

    def updateTemplate(self, templateData):
        nAtoms = self.dataDictionary['natoms']
        lineNo = 0
        while True:
            if lineNo == len(templateData):
                break
            if "LX" in templateData[lineNo]:
                templateData[lineNo] = templateData[lineNo].replace("LX", str(self.dataDictionary['lx']))
                templateData[lineNo] = templateData[lineNo].replace("LY", str(self.dataDictionary['ly']))
                templateData[lineNo] = templateData[lineNo].replace("LZ", str(self.dataDictionary['lz']))
            if "NATOMS" in templateData[lineNo]:
                templateData[lineNo] = templateData[lineNo].replace("NATOMS", str(nAtoms))
            if "NBONDS" in templateData[lineNo]:
                templateData[lineNo] = templateData[lineNo].replace("NBONDS", str(len(self.dataDictionary['bond'])))
            if "NANGLES" in templateData[lineNo]:
                templateData[lineNo] = templateData[lineNo].replace("NANGLES", str(len(self.dataDictionary['angle'])))
            if "NDIHEDRALS" in templateData[lineNo]:
                templateData[lineNo] = templateData[lineNo].replace("NDIHEDRALS", str(len(self.dataDictionary['dihedral'])))
            if "NIMPROPERS" in templateData[lineNo]:
                templateData[lineNo] = templateData[lineNo].replace("NIMPROPERS", str(len(self.dataDictionary['improper'])))
            if "<position" in templateData[lineNo]:
                if self.writeUnwrappedPositionsOnly is True:
                    for dataToWrite in list(reversed(self.dataDictionary['unwrapped_position'])):
                        stringToWrite = ''
                        for coordinate in dataToWrite:
                            stringToWrite += str(coordinate) + ' '
                        stringToWrite = stringToWrite[:-1]
                        templateData.insert(lineNo + 1, str(stringToWrite) + '\n')
                else:
                    for dataToWrite in list(reversed(self.dataDictionary['position'])):
                        stringToWrite = ''
                        for coordinate in dataToWrite:
                            stringToWrite += str(coordinate) + ' '
                        stringToWrite = stringToWrite[:-1]
                        templateData.insert(lineNo + 1, str(stringToWrite) + '\n')
            elif "<image" in templateData[lineNo]:
                if (len(self.dataDictionary['image']) != 0):
                    if self.writeUnwrappedPositionsOnly is True:
                        for atom in range(len(self.dataDictionary['image'])):
                            templateData.insert(lineNo + 1, '0 0 0\n')
                    else:
                        for dataToWrite in list(reversed(self.dataDictionary['image'])):
                            stringToWrite = ''
                            for coordinate in dataToWrite:
                                stringToWrite += str(coordinate) + ' '
                            stringToWrite = stringToWrite[:-1]
                            templateData.insert(lineNo + 1, str(stringToWrite) + '\n')
            elif "<mass" in templateData[lineNo]:
                for dataToWrite in list(reversed(self.dataDictionary['mass'])):
                    templateData.insert(lineNo + 1, str(dataToWrite) + '\n')
            elif "<diameter" in templateData[lineNo]:
                for dataToWrite in list(reversed(self.dataDictionary['diameter'])):
                    templateData.insert(lineNo + 1, str(dataToWrite) + '\n')
            elif "<type" in templateData[lineNo]:
                for dataToWrite in list(reversed(self.dataDictionary['type'])):
                    templateData.insert(lineNo + 1, str(dataToWrite) + '\n')
            elif "<body" in templateData[lineNo]:
                for dataToWrite in list(reversed(self.dataDictionary['body'])):
                    templateData.insert(lineNo + 1, str(dataToWrite) + '\n')
            elif "<bond" in templateData[lineNo]:
                for dataToWrite in list(reversed(self.dataDictionary['bond'])):
                    stringToWrite = ''
                    for bondInfo in dataToWrite:
                        stringToWrite += str(bondInfo) + ' '
                    stringToWrite = stringToWrite[:-1]
                    templateData.insert(lineNo + 1, str(stringToWrite) + '\n')
            elif "<angle" in templateData[lineNo]:
                for dataToWrite in list(reversed(self.dataDictionary['angle'])):
                    stringToWrite = ''
                    for angleInfo in dataToWrite:
                        stringToWrite += str(angleInfo) + ' '
                    stringToWrite = stringToWrite[:-1]
                    templateData.insert(lineNo + 1, str(stringToWrite) + '\n')
            elif "<dihedral" in templateData[lineNo]:
                for dataToWrite in list(reversed(self.dataDictionary['dihedral'])):
                    stringToWrite = ''
                    for dihedralInfo in dataToWrite:
                        stringToWrite += str(dihedralInfo) + ' '
                    stringToWrite = stringToWrite[:-1]
                    templateData.insert(lineNo + 1, str(stringToWrite) + '\n')
            elif "<improper" in templateData[lineNo]:
                for dataToWrite in list(reversed(self.dataDictionary['improper'])):
                    stringToWrite = ''
                    for improperInfo in dataToWrite:
                        stringToWrite += str(improperInfo) + ' '
                    stringToWrite = stringToWrite[:-1]
                    templateData.insert(lineNo + 1, str(stringToWrite) + '\n')
            elif "<charge" in templateData[lineNo]:
                for dataToWrite in list(reversed(self.dataDictionary['charge'])):
                    templateData.insert(lineNo + 1, str(dataToWrite) + '\n')
            lineNo += 1
        return templateData

    def writeData(self, newTemplateData):
        with open(self.outputFile, 'w+') as xmlFile:
            xmlFile.writelines(newTemplateData)
        print "XML file written to", self.outputFile
