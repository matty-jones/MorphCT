import copy
import numpy as np
import os
import runHoomd
import helperFunctions
import cPickle as pickle

class morphology:
    def __init__(self, morphologyXML, morphologyName, parameterDict):
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
        CGtoAAIDMaster = [] # This is a list of dictionaries. Elements in the list correspond to molecules
        # (useful for splitting out individual molecules for the xyz conversion) within the element, the
        # dictionary key is the CG site, the value is a list containing the CG type (e.g. 'thio') as the
        # first element and then another list of all the AAIDs corresponding to that CG site as the second
        # element.

        # Create a ghost particle dictionary to be added at the end of the morphology.
        # This way, we don't mess up the order of atoms later on when trying to split
        # back up into individual molecules and monomers.
        # The ghost dictionary contains all of the type T and type X particles that will
        # anchor the thiophene rings to the CG COM positions.
        ghostDictionary = {'position':[], 'image':[], 'unwrapped_position':[], 'mass':[], 'diameter':[], 'type':[], 'body':[], 'bond':[], 'angle':[], 'dihedral':[], 'improper':[], 'charge':[]}
        print "Adding molecules to the system..."
        for moleculeNumber in range(len(moleculeIDs)):
            print "Adding molecule number", moleculeNumber, "\r",
            # print "Rolling AA Index =", rollingAAIndex
            CGMoleculeDict, AAMoleculeDict, CGtoAAIDs, ghostDictionary = atomistic(moleculeIDs[moleculeNumber], self.CGDictionary, self.morphologyName, self.AATemplateFile, self.CGToTemplateAAIDs, rollingAAIndex, moleculeNumber, moleculeLengths, ghostDictionary).returnData()
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
        totalNumberOfAtoms = len(AAMorphologyDict['type']) # Should be == rollingAAIndex, but don't want to take any chances
        # Add in the wrapped positions of the ghosts. Need to know sim dims for this
        for key in ['lx', 'ly', 'lz']:
            ghostDictionary[key] = AAMorphologyDict[key]
        ghostDictionary = helperFunctions.addWrappedPositions(ghostDictionary)
        for key in ['lx', 'ly', 'lz']:
            ghostDictionary.pop(key)
#        # THIS HAS NOW ALREADY BEEN DONE WHEN THE ATOMIDs WERE INCREMENTED PREVIOUSLY
#        # Increment the bond IDs to match the ghost particle IDs
#        for bondNo, bond in enumerate(ghostDictionary['bond']):
#            if ('C4' not in bond[0]) and ('C7' not in bond[0]):
#                ghostDictionary['bond'][bondNo] = [bond[0], bond[1]+totalNumberOfAtoms, bond[2]+totalNumberOfAtoms]
#            else:
#                ghostDictionary['bond'][bondNo] = [bond[0], bond[1]+totalNumberOfAtoms, bond[2]]
        # Now append all ghosts to morphology
        for key in ghostDictionary.keys():
            AAMorphologyDict[key] += ghostDictionary[key]
        # Finally, update the number of atoms
        AAMorphologyDict['natoms'] += len(ghostDictionary['type'])
        print "\n"
        print "Writing XML file..."
        AAFileName = './outputFiles/'+morphologyName+'/morphology/'+morphologyName+'.xml'
        writeXML(AAMorphologyDict, './templates/template.xml', AAFileName)
        toPickle = (AAFileName, CGMorphologyDict, AAMorphologyDict, CGtoAAIDMaster, [], boxSize)
        print "Writing pickle file..."
        pickleFileName = './outputFiles/'+morphologyName+'/morphology/'+morphologyName+'.pickle'
        with open(pickleFileName, 'w+') as pickleFile:
            pickle.dump(toPickle, pickleFile)
        print "Pickle file written to", pickleFileName
        return AAFileName, CGMorphologyDict, AAMorphologyDict, CGtoAAIDMaster, [], boxSize
        # exit()
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
        #morphologySegments = []
        #for moleculeAtoms in moleculeIDs:
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
            while addedNewSite == True:
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
            moleculeLength = np.round(currentNumberOfSitesInMolecule/len(self.CGToTemplateAAIDs.keys())) # Just in case of floating point errors
            moleculeLengths.append(moleculeLength)
        return moleculeList, moleculeLengths


class molecule:
    def __init__(self, moleculeIDs, CGDictionary):
        #### TOLERANCE SHOULD BE IN A PARAMETER FILE ####
        self.tolerance = np.pi/6.
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
            if firstAtomInSegment == True:
                axisVector = helperFunctions.findAxis(atomUnderConsiderationCoords, neighbouringBackboneAtomCoords)
                firstAtomInSegment = False
                print "No longer first atom in segment"
            separationVector = helperFunctions.findAxis(atomUnderConsiderationCoords, neighbouringBackboneAtomCoords)
            dotProduct = np.dot(axisVector, separationVector)
            print dotProduct
            #print "Axis Vector for this atom =", axisVector, "Separation Vector for this atom =", separationVector, "Dot Product for this atom =", dotProduct
            if abs(dotProduct-1.0) <= 1E-8:
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
    def __init__(self, siteIDs, CGDictionary, morphologyName, templateFile, CGToTemplateAAIDs, rollingAAIndex, moleculeIndex, moleculeLengths, ghostDictionary):
        # This class sees individual molecules.
        self.noAtomsInMorphology = rollingAAIndex
        self.moleculeIndex = moleculeIndex
        self.moleculeLengths = moleculeLengths
        self.siteIDs = siteIDs
        self.CGDictionary = CGDictionary
        self.CGMonomerDictionary = self.getCGMonomerDict()
        AATemplateDictionary = helperFunctions.loadMorphologyXML(templateFile)
        self.AATemplateDictionary = helperFunctions.addUnwrappedPositions(AATemplateDictionary)
        # thioAA, alk1AA, alk2AA are now all unused. Instead, use the CGToTemplateAAIDs dictionary which can have arbitrary length
        self.CGToTemplateAAIDs = CGToTemplateAAIDs
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
        CGMonomerDictionary = {'position':[], 'image':[], 'mass':[], 'diameter':[], 'type':[], 'body':[], 'bond':[], 'angle':[], 'dihedral':[], 'improper':[], 'charge':[], 'lx':0, 'ly':0, 'lz':0}
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
        AADictionary = {'position':[], 'image':[], 'unwrapped_position':[], 'mass':[], 'diameter':[], 'type':[], 'body':[], 'bond':[], 'angle':[], 'dihedral':[], 'improper':[], 'charge':[], 'lx':0, 'ly':0, 'lz':0}
        CGCoMs = self.getAATemplatePosition(self.CGToTemplateAAIDs)
        # Need to keep track of the atom ID numbers globally - runFineGrainer sees individual monomers, atomistic sees molecules and the XML needs to contain the entire morphology.
        noAtomsInMolecule = 0
        CGTypeList = {}
        for siteID in self.siteIDs:
            CGTypeList[siteID] = self.CGDictionary['type'][siteID]
        monomerList = self.sortIntoMonomers(CGTypeList)
        currentMonomerIndex = sum(self.moleculeLengths[:self.moleculeIndex])
        atomIDLookupTable = {}
        totalPermittedAtoms = len(self.AATemplateDictionary['type'])*len(monomerList)
        for monomerNo, monomer in enumerate(monomerList):
            thisMonomerDictionary = copy.deepcopy(self.AATemplateDictionary)
            for key in ['lx', 'ly', 'lz']:
                thisMonomerDictionary[key] = self.CGDictionary[key]
            if len(thisMonomerDictionary['image']) == 0:
                thisMonomerDictionary['image'] = [[0, 0, 0]]*len(thisMonomerDictionary['position'])
            for siteID in monomer:
                sitePosn = np.array(self.CGDictionary['unwrapped_position'][siteID])
                siteTranslation = sitePosn - CGCoMs[CGTypeList[siteID]]
                # See if it works without rotating the groups
                ## Now rotate the functional groups in place so that they connect more nicely
                #thioAlk1Axis = helperFunctions.findAxis(thioPosn, alk1Posn)
                #alk1Alk2Axis = helperFunctions.findAxis(alk1Posn, alk2Posn)
                #thisMonomerDictionary = self.rotateFunctionalGroups(thisMonomerDictionary, thioAlk1Axis, alk1Alk2Axis, thioAACOM, alk1AACOM, alk2AACOM, thioAAIDs, alk1AAIDs, alk2AAIDs)
                atomIDLookupTable[siteID] = [CGTypeList[siteID], [x + noAtomsInMolecule + self.noAtomsInMorphology for x in self.CGToTemplateAAIDs[CGTypeList[siteID]]]]
                # Add the atoms in based on the CG site position
                for AAID in self.CGToTemplateAAIDs[CGTypeList[siteID]]:
                    thisMonomerDictionary['unwrapped_position'][AAID] = list(np.array(thisMonomerDictionary['unwrapped_position'][AAID])+siteTranslation)
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
                    ghostDictionary['unwrapped_position'].append(helperFunctions.calcCOM(AAIDPositions, listOfAtomTypes = AAIDAtomTypes))
                    ghostDictionary['mass'].append(1.0)
                    ghostDictionary['diameter'].append(1.0)
                    ghostDictionary['type'].append('R'+str(CGTypeList[siteID]))
                    ghostDictionary['body'].append(currentMonomerIndex)
                    ghostDictionary['charge'].append(0.0)
                    # Then create the corresponding CG anchorpoint
                    ghostDictionary['unwrapped_position'].append(self.CGDictionary['unwrapped_position'][siteID])
                    ghostDictionary['mass'].append(1.0)
                    ghostDictionary['diameter'].append(1.0)
                    ghostDictionary['type'].append('X'+str(CGTypeList[siteID]))
                    ghostDictionary['body'].append(-1)
                    ghostDictionary['charge'].append(0.0)
                    # Now create a bond between them
                    # We want to bond together the previous two ghost particles, so this should work as it requires no knowledge of the number of ghost particles already in the system.
                    ghostDictionary['bond'].append([str(ghostDictionary['type'][-2])+'-'+str(ghostDictionary['type'][-1]), len(ghostDictionary['type'])-2, len(ghostDictionary['type'])-1])
                else:
                    # Create a ghost particle that describe the CG anchorpoint for the non-rigid body
                    ghostDictionary['unwrapped_position'].append(self.CGDictionary['unwrapped_position'][siteID])
                    ghostDictionary['mass'].append(1.0)
                    ghostDictionary['diameter'].append(1.0)
                    ghostDictionary['type'].append('X'+str(CGTypeList[siteID]))
                    ghostDictionary['body'].append(-1)
                    ghostDictionary['charge'].append(0.0)
                    # Add in bonds between the CG anchorpoints and the atom closest to the CG site
                    # Find the atom closest to the CG site
                    closestAtomID = None
                    closestAtomPosn = 1E99
                    for AAID, AAPosition in enumerate(thisMonomerDictionary['unwrapped_positions']):
                        separation = helperFunctions.calculateSeparation(self.CGDictionary['unwrapped_position'][siteID], AAPosition)
                        if separation < closestAtomPosn:
                            closestAtomPosn = separation
                            closestAtomID = AAID
                    # Add in the bond:
                    # Note that, in order to distinguish between the ghostAtomIDs and the realAtomIDs, I've put an underscore in front of the closestAtomID. When adding the ghost dictionary to to the system later, this will be unchanged. However, these values will need to be incremented when incrementing the realAtomIDs for the monomer and the molecule.
                    ghostDictionary['bond'].append([str(ghostDictionary['type'][-1])+'-'+str(thisMonomerDictionary['type'][closestAtomID]), len(ghostDictionary['type'])-1, '_'+str(closestAtomID)])
                    # I'm not sure how this is going to work later when I incorporate the ghost dictionary and all of the atom ID's change. Perhaps I can change this one to '_'+closestAtomID and put in a check later not to change the ID's that have an underscore?
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
                if atFinalMonomer == True:
                    break
                # Work out which key to write the constraint to based on its length:
                # 3 = Bond, 4 = Angle, 5 = Dihedral, 6 = Improper
                constraintType = ['bond', 'angle', 'dihedral', 'improper']
                thisMonomerDictionary[constraintType[len(constraint)-3]].append(constraint)
            # Finally, increment the atom IDs and then update the AADictionary.
            thisMonomerDictionary, ghostDictionary = helperFunctions.incrementAtomIDs(thisMonomerDictionary, ghostDictionary, noAtomsInMolecule)
            # Find the connecting atoms to the terminating units based on monomer number
            if len(self.moleculeTerminatingConnections) != 0:
                if monomerNo == 0:
                    startAtomIndex = noAtomsInMolecule + self.moleculeTerminatingConnections[0][1]
                elif monomerNo == len(monomerList)-1:
                    endAtomIndex = noAtomsInMolecule + self.moleculeTerminatingConnections[1][1]
            noAtomsInMolecule += 25
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
        AADictionary, ghostDictionary = helperFunctions.incrementAtomIDs(AADictionary, ghostDictionary, self.noAtomsInMorphology)
        return AADictionary, atomIDLookupTable, ghostDictionary


#            atomIDLookupTable[thioID] = ['thio', [x + noAtomsInMolecule + self.noAtomsInMorphology for x in thioAAIDs]]
#            atomIDLookupTable[alk1ID] = ['alk1', [x + noAtomsInMolecule + self.noAtomsInMorphology for x in alk1AAIDs]]
#            atomIDLookupTable[alk2ID] = ['alk2', [x + noAtomsInMolecule + self.noAtomsInMorphology for x in alk2AAIDs]]
#            for thioAAID in thioAAIDs:
#                thisMonomerDictionary['unwrapped_position'][thioAAID] = list(np.array(thisMonomerDictionary['unwrapped_position'][thioAAID])+thioTranslation)
#                # Set the thiophenes to be rigid bodies
#                if thisMonomerDictionary['type'][thioAAID] == 'H1':
#                    thisMonomerDictionary['body'][thioAAID] = -1
#                else:
#                    thisMonomerDictionary['body'][thioAAID] = currentMonomerIndex
#            for alk1AAID in alk1AAIDs:
#                thisMonomerDictionary['unwrapped_position'][alk1AAID] = list(np.array(thisMonomerDictionary['unwrapped_position'][alk1AAID])+alk1Translation)
#            for alk2AAID in alk2AAIDs:
#                thisMonomerDictionary['unwrapped_position'][alk2AAID] = list(np.array(thisMonomerDictionary['unwrapped_position'][alk2AAID])+alk2Translation)
#            thisMonomerDictionary = helperFunctions.addWrappedPositions(thisMonomerDictionary)
#            thisMonomerDictionary = helperFunctions.addMasses(thisMonomerDictionary)
#            thisMonomerDictionary = helperFunctions.addDiameters(thisMonomerDictionary)
#            
#            
#            #### DEBUG CODE ####
#            # listOfMasses = []
#            # listOfPositions = []
#            # for atomNo in thioAAIDs:
#            #     listOfPositions.append(thisMonomerDictionary['unwrapped_position'][atomNo])
#            #     if ('H' in thisMonomerDictionary['type'][atomNo]):
#            #         listOfMasses.append(1.0)
#            #     elif ('C' in thisMonomerDictionary['type'][atomNo]):
#            #         listOfMasses.append(12.0)
#            #     elif ('S' in thisMonomerDictionary['type'][atomNo]):
#            #         listOfMasses.append(32.0)
#            # print len(listOfPositions), len(listOfMasses)
#            # COMPosn = calcCOM(listOfPositions, listOfMasses)
#            # print "CGCOM =", thioPosn
#            # print "AACOM =", COMPosn
#            #####################
#
#
#            # Now add in the two ghost particles for the thiophene ring:
#            # Type T = COM ghost particle of the thiophene ring that has the same rigid body as the thiophene
#            #          atoms, but does not interact with them. Its initial position is has the COM from the
#            #          CG morph.
#            # Type X1 = Anchor particle that is never integrated over, has no interactions with anything except
#            #          a bond to the corresponding type T. The type X particle never moves from its position,
#            #          which is the Type A position from the CG morph (COM of Thiophene ring)
#
#            # Type T:
#            ghostDictionary['unwrapped_position'].append(list(thioPosn))
#            ghostDictionary['mass'].append(1.0)
#            ghostDictionary['diameter'].append(1.0)
#            ghostDictionary['type'].append('T')
#            ghostDictionary['body'].append(currentMonomerIndex)
#            ghostDictionary['charge'].append(0.0)
#
#            # Type X1:
#            ghostDictionary['unwrapped_position'].append(list(thioPosn))
#            ghostDictionary['mass'].append(1.0)
#            ghostDictionary['diameter'].append(1.0)
#            ghostDictionary['type'].append('X1')
#            ghostDictionary['body'].append(-1)
#            ghostDictionary['charge'].append(0.0)
#
#
#            # Now add the T-X1 bond
#            ghostDictionary['bond'].append(['T-X1', 4*currentMonomerIndex, 4*currentMonomerIndex+1])
#            # NOTE: The 4x is because I am adding 4 different ghost particles to the system
#            # It keeps the numbering here constant. Later on it gets added on
#            
#
#            # Now add in the two ghost particles for the alkyl sidechains:
#            # Type X2 = Anchor particle that is never integrated over, has no interactions with anything except
#            #           except a bond to the alkyl sidechain atom C4 to constrain the alk1 COM posn.
#            #           Its initial position is has the COM from the CG morph.
#            # Type X3 = Anchor particle that is never integrated over, has no interactions with anything except
#            #           except a bond to the alkyl sidechain atom C7 to constrain the alk2 COM posn.
#            #           Its initial position is has the COM from the CG morph.
#
#            # Type X2:
#            ghostDictionary['unwrapped_position'].append(list(alk1Posn))
#            ghostDictionary['mass'].append(1.0)
#            ghostDictionary['diameter'].append(1.0)
#            ghostDictionary['type'].append('X2')
#            ghostDictionary['body'].append(-1)
#            ghostDictionary['charge'].append(0.0)
#
#            for atomID, atomType in enumerate(thisMonomerDictionary['type']):
#                if atomType == 'C4':
#                    C4Index = atomID
#                if atomType == 'C7':
#                    C7Index = atomID
#
#            # Now add the X2-C4 bond
#            ghostDictionary['bond'].append(['X2-C4', 4*currentMonomerIndex+2, self.noAtomsInMorphology + noAtomsInMolecule + C4Index])
#
#            # Type X3:
#            ghostDictionary['unwrapped_position'].append(list(alk2Posn))
#            ghostDictionary['mass'].append(1.0)
#            ghostDictionary['diameter'].append(1.0)
#            ghostDictionary['type'].append('X3')
#            ghostDictionary['body'].append(-1)
#            ghostDictionary['charge'].append(0.0)
#
#            # Now add the X3-C7 bond
#            ghostDictionary['bond'].append(['X3-C7', 4*currentMonomerIndex+3, self.noAtomsInMorphology + noAtomsInMolecule + C7Index])
#
#
#            # Now add in the inter-monomer bond.
#            # In the template, the inter-monomer bonding carbons are given by type C10 and C1, which are thioAAIDs 0 and 3
#            if thioID != moleculeEnds[1]:
#                # We are at the start of the molecule, or somewhere in the middle and so an intermonomer bond is required
#                # to connect to the next monomer. This bond is between thioAAID 3 and thioAAID 0+25 (again the 25 is hard
#                # coded for this template, because it has 25 atoms in it)
#                # Create the bond (1 bond)
#                thisMonomerDictionary['bond'].append(['C1-C10', 3, 25])
#                # Create the angles (4 angles)
#                thisMonomerDictionary['angle'].append(['C1-C10-C9', 3, 25, 26])
#                thisMonomerDictionary['angle'].append(['C1-C10-S1', 3, 25, 29])
#                thisMonomerDictionary['angle'].append(['S1-C1-C10', 4, 3, 25])
#                thisMonomerDictionary['angle'].append(['C2-C1-C10', 2, 3, 25])                
#                # Create the dihedrals that are necessary to define the forcefield (the commented out ones are not required)
#                # thisMonomerDictionary['dihedral'].append(['C1-C10-C9-H1', 3, 25, 26, 49])
#                thisMonomerDictionary['dihedral'].append(['C1-C10-C9-C2', 3, 25, 26, 27])
#                thisMonomerDictionary['dihedral'].append(['C1-C10-S1-C1', 3, 25, 29, 28])
#                # thisMonomerDictionary['dihedral'].append(['C10-C1-S1-C1', 25, 3, 4, 0])
#                # thisMonomerDictionary['dihedral'].append(['C10-C1-C2-C9', 25, 3, 2, 1])
#                # thisMonomerDictionary['dihedral'].append(['C10-C1-C2-C3', 25, 3, 2, 5])
#                # thisMonomerDictionary['dihedral'].append(['S1-C1-C10-C9', 4, 3, 25, 26])
#                thisMonomerDictionary['dihedral'].append(['S1-C1-C10-S1', 4, 3, 25, 29])
#                # thisMonomerDictionary['dihedral'].append(['C2-C1-C10-C9', 2, 3, 25, 26])
#                thisMonomerDictionary['dihedral'].append(['C2-C1-C10-S1', 2, 3, 25, 29])
#                
#            # Now we update the atom IDs to mirror the fact that we have added an additional monomer to the system
#            thisMonomerDictionary = helperFunctions.incrementAtomIDs(thisMonomerDictionary, noAtomsInMolecule)
#            noAtomsInMolecule += 25
#            currentMonomerIndex += 1
#            AADictionary = self.updateMoleculeDictionary(thisMonomerDictionary, AADictionary)
#        AADictionary['natoms'] = noAtomsInMolecule
#        for key in ['lx', 'ly', 'lz']:
#            AADictionary[key] = thisMonomerDictionary[key]
#        AADictionary, terminatingHydrogenIDs = helperFunctions.addTerminatingHydrogens(AADictionary)
#        # Update the atomID CG LookupTable with the terminating Hydrogens
#        atomIDLookupTable[moleculeEnds[0]][1].append(terminatingHydrogenIDs[0] + self.noAtomsInMorphology)
#        atomIDLookupTable[moleculeEnds[1]][1].append(terminatingHydrogenIDs[1] + self.noAtomsInMorphology)
#        # Now the molecule is done, we need to add on the correct identifying numbers for all the bonds, angles and dihedrals
#        # (just as we did between monomers) for the other molecules in the system, so that they all connect to the right atoms
#        AADictionary = helperFunctions.incrementAtomIDs(AADictionary, self.noAtomsInMorphology)
#        return AADictionary, atomIDLookupTable, ghostDictionary

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
            while addedNewSite == True:
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


    def rotateFunctionalGroups(self, atomDictionary, thioAlk1Axis, alk1Alk2Axis, thioAACOM, alk1AACOM, alk2AACOM, thioAAIDs, alk1AAIDs, alk2AAIDs):
        # Treat Thiophene by rotating it so that the COM-C2 axis is the same as the thio-alk1 axis
        # Determine the rotation we want
        C2AtomID = atomDictionary['type'].index('C2')
        COM_C2Axis = helperFunctions.findAxis(thioAACOM, atomDictionary['unwrapped_position'][C2AtomID])
        thioRotation = helperFunctions.rotationMatrix(COM_C2Axis, thioAlk1Axis)
        # Now apply rotation
        for thiopheneAtomID in thioAAIDs:
            # First, move all thiophene atoms to the origin for rotation
            newPosn = np.array(atomDictionary['unwrapped_position'][thiopheneAtomID]) - np.array(thioAACOM)
            # Do the matrix multiplication for the rotation
            newPosn = np.transpose(thioRotation*np.transpose(np.matrix(newPosn)))
            # Put the atom back in the correct COM position
            newPosn += np.array(thioAACOM)
            # Update the atom dictionary
            atomDictionary['unwrapped_position'][thiopheneAtomID] = [newPosn[0,0], newPosn[0,1], newPosn[0,2]]
        # Treat Alk1 by rotating it so that the C5-C3 axis is the same as the alk1-thio axis
        C3AtomID = atomDictionary['type'].index('C3')
        C5AtomID = atomDictionary['type'].index('C5')
        C5_C3Axis = helperFunctions.findAxis(atomDictionary['unwrapped_position'][C5AtomID], atomDictionary['unwrapped_position'][C3AtomID])
        alk1Rotation = helperFunctions.rotationMatrix(C5_C3Axis, -thioAlk1Axis)
        # Now apply rotation
        for alk1AtomID in alk1AAIDs:
            newPosn = np.array(atomDictionary['unwrapped_position'][alk1AtomID]) - np.array(alk1AACOM)
            newPosn = np.transpose(alk1Rotation*np.transpose(np.matrix(newPosn)))
            newPosn += np.array(alk1AACOM)
            atomDictionary['unwrapped_position'][alk1AtomID] = [newPosn[0,0], newPosn[0,1], newPosn[0,2]]
        # Treat Alk2 by rotating it so that the C8-C6 axis is the same as the alk2-alk1 axis
        C6AtomID = atomDictionary['type'].index('C6')
        C8AtomID = atomDictionary['type'].index('C8')
        C8_C6Axis = helperFunctions.findAxis(atomDictionary['unwrapped_position'][C8AtomID], atomDictionary['unwrapped_position'][C6AtomID])
        alk2Rotation = helperFunctions.rotationMatrix(C8_C6Axis, -alk1Alk2Axis)
        # Now apply rotation
        for alk2AtomID in alk2AAIDs:
            newPosn = np.array(atomDictionary['unwrapped_position'][alk2AtomID]) - np.array(alk2AACOM)
            newPosn = np.transpose(alk2Rotation*np.transpose(np.matrix(newPosn)))
            newPosn += np.array(alk2AACOM)
            atomDictionary['unwrapped_position'][alk2AtomID] = [newPosn[0,0], newPosn[0,1], newPosn[0,2]]
        return atomDictionary


    def identifyMonomerSites(self, backboneID):
        # 3 length list of the form: [thio, bonded-alk1, bonded-alk2]
        monomerSites = [backboneID]
        for bond in self.CGDictionary['bond']:
            if (bond[0] == 'bondB') and (bond[1] == backboneID):
                monomerSites.append(bond[2])
                break
            elif (bond[0] == 'bondB') and (bond[2] == backboneID):
                monomerSites.append(bond[1])
                break
        for bond in self.CGDictionary['bond']:
            if (bond[0] == 'bondC') and (bond[1] == monomerSites[-1]):
                monomerSites.append(bond[2])
                break
            elif (bond[0] == 'bondC') and (bond[2] == monomerSites[-1]):
                monomerSites.append(bond[1])
                break
        return monomerSites


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
            CGCoMs[siteName] = helperFunctions.calcCOM(sitePositions, listOfAtomTypes = siteTypes)
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
                if self.writeUnwrappedPositionsOnly == True:
                    for dataToWrite in list(reversed(self.dataDictionary['unwrapped_position'])):
                        stringToWrite = ''
                        for coordinate in dataToWrite:
                            stringToWrite += str(coordinate)+' '
                        stringToWrite = stringToWrite[:-1]
                        templateData.insert(lineNo+1, str(stringToWrite)+'\n')
                else:
                    for dataToWrite in list(reversed(self.dataDictionary['position'])):
                        stringToWrite = ''
                        for coordinate in dataToWrite:
                            stringToWrite += str(coordinate)+' '
                        stringToWrite = stringToWrite[:-1]
                        templateData.insert(lineNo+1, str(stringToWrite)+'\n')
            elif "<image" in templateData[lineNo]:
                if (len(self.dataDictionary['image']) != 0):
                    if self.writeUnwrappedPositionsOnly == True:
                        for atom in range(len(self.dataDictionary['image'])):
                            templateData.insert(lineNo+1, '0 0 0\n')
                    else:
                        for dataToWrite in list(reversed(self.dataDictionary['image'])):
                            stringToWrite = ''
                            for coordinate in dataToWrite:
                                stringToWrite += str(coordinate)+' '
                            stringToWrite = stringToWrite[:-1]
                            templateData.insert(lineNo+1, str(stringToWrite)+'\n')
            elif "<mass" in templateData[lineNo]:
                for dataToWrite in list(reversed(self.dataDictionary['mass'])):
                    templateData.insert(lineNo+1, str(dataToWrite)+'\n')
            elif "<diameter" in templateData[lineNo]:
                for dataToWrite in list(reversed(self.dataDictionary['diameter'])):
                    templateData.insert(lineNo+1, str(dataToWrite)+'\n')
            elif "<type" in templateData[lineNo]:
                for dataToWrite in list(reversed(self.dataDictionary['type'])):
                    templateData.insert(lineNo+1, str(dataToWrite)+'\n')
            elif "<body" in templateData[lineNo]:
                for dataToWrite in list(reversed(self.dataDictionary['body'])):
                    templateData.insert(lineNo+1, str(dataToWrite)+'\n')
            elif "<bond" in templateData[lineNo]:
                for dataToWrite in list(reversed(self.dataDictionary['bond'])):
                    stringToWrite = ''
                    for bondInfo in dataToWrite:
                        stringToWrite += str(bondInfo)+' '
                    stringToWrite = stringToWrite[:-1]
                    templateData.insert(lineNo+1, str(stringToWrite)+'\n')
            elif "<angle" in templateData[lineNo]:
                for dataToWrite in list(reversed(self.dataDictionary['angle'])):
                    stringToWrite = ''
                    for angleInfo in dataToWrite:
                        stringToWrite += str(angleInfo)+' '
                    stringToWrite = stringToWrite[:-1]
                    templateData.insert(lineNo+1, str(stringToWrite)+'\n')
            elif "<dihedral" in templateData[lineNo]:
                for dataToWrite in list(reversed(self.dataDictionary['dihedral'])):
                    stringToWrite = ''
                    for dihedralInfo in dataToWrite:
                        stringToWrite += str(dihedralInfo)+' '
                    stringToWrite = stringToWrite[:-1]
                    templateData.insert(lineNo+1, str(stringToWrite)+'\n')
            elif "<improper" in templateData[lineNo]:
                for dataToWrite in list(reversed(self.dataDictionary['improper'])):
                    stringToWrite = ''
                    for improperInfo in dataToWrite:
                        stringToWrite += str(improperInfo)+' '
                    stringToWrite = stringToWrite[:-1]
                    templateData.insert(lineNo+1, str(stringToWrite)+'\n')
            elif "<charge" in templateData[lineNo]:
                for dataToWrite in list(reversed(self.dataDictionary['charge'])):
                    templateData.insert(lineNo+1, str(dataToWrite)+'\n')
            lineNo += 1
        return templateData
            
    def writeData(self, newTemplateData):
        with open(self.outputFile, 'w+') as xmlFile:
            xmlFile.writelines(newTemplateData)
        print "XML file written to", self.outputFile
