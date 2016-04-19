import copy
import numpy as np
import os
import runHoomd
import helperFunctions
import cPickle as pickle

class morphology:
    def __init__(self, morphologyName, sigma=1.):
        # Sigma is the `compression value' in Angstroms that has been used to scale the morphology
        # E.G. the P3HT template uses sigma = 1, but the Marsh morphologies use sigma = 3.
        self.xmlPath = str(morphologyName)
        self.CGDictionary = helperFunctions.loadMorphologyXML(self.xmlPath, sigma=sigma)
        self.CGDictionary = helperFunctions.addUnwrappedPositions(self.CGDictionary)


    def analyseMorphology(self):
        print "Finding molecules..."
        moleculeIDs = self.splitMolecules()
        # Obtain the morphology name:
        slashLocs = helperFunctions.findIndex(self.xmlPath, '/')
        morphologyName = self.xmlPath[slashLocs[-1]+1:-4]
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
        print "Adding molecules to the system..."
        for moleculeNumber in range(len(moleculeIDs)):
            print "Adding molecule number", moleculeNumber, "\r",
            # print "Rolling AA Index =", rollingAAIndex
            CGMoleculeDict, AAMoleculeDict, CGtoAAIDs = atomistic(moleculeIDs[moleculeNumber], self.CGDictionary, morphologyName, rollingAAIndex).returnData()
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
        bondList = copy.deepcopy(self.CGDictionary['bond'])
        while len(bondList) > 0:
            # Add the first two atoms of the molecule
            thisMolecule = [bondList[0][1], bondList[0][2]]
            addedNewAtom = True
            while addedNewAtom == True:
                addedNewAtom = False
                bondPopList = []
                for bond in bondList:
                    if (bond[1] in thisMolecule) and (bond[2] not in thisMolecule):
                        thisMolecule.append(bond[2])
                        addedNewAtom = True
                    elif (bond[2] in thisMolecule) and (bond[1] not in thisMolecule):
                        thisMolecule.append(bond[1])
                        addedNewAtom = True
                    elif (bond[1] in thisMolecule) and (bond[2] in thisMolecule):
                        pass
                    else:
                        continue
                    bondPopList.append(bondList.index(bond))
                bondPopList.sort(reverse=True)
                for bondIndex in bondPopList:
                    bondList.pop(bondIndex)
            moleculeList.append(thisMolecule)
        return moleculeList


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
    def __init__(self, atomIDs, CGDictionary, morphologyName, rollingAAIndex):
        self.noAtomsInMorphology = rollingAAIndex
        self.atomIDs = atomIDs
        self.CGDictionary = CGDictionary
        self.CGMonomerDictionary = self.getCGMonomerDict()
        polymerBackboneIDs, moleculeEnds = self.obtainBackboneAtoms()
        # THIS BIT IS SPECIFIC TO P3HT SO ADDITIONAL MODELS WILL NEED TO BE
        # HARD CODED IN HERE
        templateDict, thioAA, alk1AA, alk2AA = self.loadAATemplate()
        self.AATemplateDictionary = templateDict
        self.AADictionary, self.atomIDLookupTable = self.runFineGrainer(thioAA, alk1AA, alk2AA, polymerBackboneIDs, moleculeEnds)

    def returnData(self):
        return self.CGMonomerDictionary, self.AADictionary, self.atomIDLookupTable

                
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

    
    def loadAATemplate(self):
        templateDict = morphology('./templates/mid3HT.xml').CGDictionary
        templateDict = self.trimTemplateDihedrals(templateDict)
        # The thiophene ring consists of S1, C10, C9, C2, C1, H1 bonded to C9
        thioAtomIDs = []
        alk1AtomIDs = []
        alk2AtomIDs = []
        # Only one sulfur (S1) and alkyl bonding carbon (C2) so we can use index here
        thioAtomIDs.append(templateDict['type'].index('S1'))
        thioAtomIDs.append(templateDict['type'].index('C2'))
        thioAtomIDs.append(templateDict['type'].index('C1'))
        thioAtomIDs.append(templateDict['type'].index('C10'))
        # By examining the bonds, we can add all of the correct atoms to the correct lists.
        for bond in templateDict['bond']:
            if ('C9' in bond[0]) and ('H1' in bond[0]):
                if bond[1] not in thioAtomIDs:
                    thioAtomIDs.append(bond[1])
                if bond[2] not in thioAtomIDs:
                    thioAtomIDs.append(bond[2])
            elif ('H1' in bond[0]) and (('C3' in bond[0]) or ('C4' in bond[0]) or ('C5' in bond[0])):
                if bond[1] not in alk1AtomIDs:
                    alk1AtomIDs.append(bond[1])
                if bond[2] not in alk1AtomIDs:
                    alk1AtomIDs.append(bond[2])
            elif ('H1' in bond[0]) and (('C6' in bond[0]) or ('C7' in bond[0]) or ('C8' in bond[0])):
                if bond[1] not in alk2AtomIDs:
                    alk2AtomIDs.append(bond[1])
                if bond[2] not in alk2AtomIDs:
                    alk2AtomIDs.append(bond[2])
        thioAtomIDs.sort()
        alk1AtomIDs.sort()
        alk2AtomIDs.sort()
        return templateDict, thioAtomIDs, alk1AtomIDs, alk2AtomIDs


    def trimTemplateDihedrals(self, templateDict):
        # Remove the dihedrals that are not required in order to properly define the P3HT forcefield
        # (hardcoded to P3HT)
        popList = []
        for dihedralNo in range(len(templateDict['dihedral'])):
            dihedralAtoms = templateDict['dihedral'][dihedralNo][0]
            # Don't need any dihedrals including H
            if 'H1' in dihedralAtoms:
                popList.append(dihedralNo)
            # Some other dihedrals are overdefined. Remove them here:
            if (dihedralAtoms == 'C9-C2-C1-S1') or (dihedralAtoms == 'C1-C2-C3-C4') or (dihedralAtoms == 'C3-C2-C1-S1'):
                popList.append(dihedralNo)
        popList.sort(reverse=True)
        for popElement in popList:
            templateDict['dihedral'].pop(popElement)
        return templateDict

    
    def getCGMonomerDict(self):
        CGMonomerDictionary = {'position':[], 'image':[], 'velocity':[], 'mass':[], 'diameter':[], 'type':[], 'body':[], 'bond':[], 'angle':[], 'dihedral':[], 'improper':[], 'charge':[], 'lx':0, 'ly':0, 'lz':0}
        # First, do just the positions and find the newAtomIDs for each CG site
        for atomID in self.atomIDs:
            CGMonomerDictionary['position'].append(self.CGDictionary['position'][atomID])
        # Now sort out the other one-per-atom properties
        for key in ['image', 'velocity', 'mass', 'diameter', 'type', 'body', 'charge']:
            if len(self.CGDictionary[key]) != 0:
                for atomID in self.atomIDs:
                    CGMonomerDictionary[key].append(self.CGDictionary[key][atomID])
        # Now rewrite the bonds based on the newAtomIDs
        for key in ['bond', 'angle', 'dihedral', 'improper']:
            for element in self.CGDictionary[key]:
                for atomID in self.atomIDs:
                    if (atomID in element) and (element not in CGMonomerDictionary[key]):
                        CGMonomerDictionary[key].append(element)
        # Now update the box parameters
        for key in ['lx', 'ly', 'lz']:
            CGMonomerDictionary[key] = self.CGDictionary[key]
        CGMonomerDictionary = helperFunctions.addUnwrappedPositions(CGMonomerDictionary)
        CGMonomerDictionary['natoms'] = len(CGMonomerDictionary['position'])
        return CGMonomerDictionary

    
    def runFineGrainer(self, thioAAIDs, alk1AAIDs, alk2AAIDs, polymerBackboneIDs, moleculeEnds):
        AADictionary = {'position':[], 'image':[], 'unwrapped_position':[], 'velocity':[], 'mass':[], 'diameter':[], 'type':[], 'body':[], 'bond':[], 'angle':[], 'dihedral':[], 'improper':[], 'charge':[], 'lx':0, 'ly':0, 'lz':0}
        thioAACOM, alk1AACOM, alk2AACOM = self.getAATemplatePosition(thioAAIDs, alk1AAIDs, alk2AAIDs)
        atomIDLookupTable = {}
        # Begin at one terminating monomer and build up the monomers
        bondedCGSites = []
        for backboneID in polymerBackboneIDs:
            monomerSites = self.identifyMonomerSites(backboneID)
            bondedCGSites.append(monomerSites)
        noAtomsInMolecule = 0
        # Normalisation no longer needed, but need to keep track of the atom ID numbers globally - runFineGrainer sees individual monomers, atomistic sees molecules and the XML needs to contain the entire morphology.
        # If this isn't the first molecule in the morphology, we need to offset the CG indices in the atomIDLookupTable, because runhoomd treats each molecule as isolated
        #CGIDOffset = bondedCGSites[0][0]
        for monomer in bondedCGSites:
            thisMonomerDictionary = copy.deepcopy(self.AATemplateDictionary)
            for key in ['lx', 'ly', 'lz']:
                thisMonomerDictionary[key] = self.CGDictionary[key]
            if len(thisMonomerDictionary['velocity']) == 0:
                thisMonomerDictionary['velocity'] = [[0.0, 0.0, 0.0]]*len(thisMonomerDictionary['position'])
            if len(thisMonomerDictionary['image']) == 0:
                thisMonomerDictionary['image'] = [[0, 0, 0]]*len(thisMonomerDictionary['position'])
            thioID = monomer[0]
            alk1ID = monomer[1]
            alk2ID = monomer[2]
            thioPosn = np.array(self.CGDictionary['unwrapped_position'][thioID])
            alk1Posn = np.array(self.CGDictionary['unwrapped_position'][alk1ID])
            alk2Posn = np.array(self.CGDictionary['unwrapped_position'][alk2ID])
            thioTranslation = thioPosn - thioAACOM
            alk1Translation = alk1Posn - alk1AACOM
            alk2Translation = alk2Posn - alk2AACOM
            # Now rotate the functional groups in place so that they connect more nicely
            thioAlk1Axis = helperFunctions.findAxis(thioPosn, alk1Posn)
            alk1Alk2Axis = helperFunctions.findAxis(alk1Posn, alk2Posn)
            thisMonomerDictionary = self.rotateFunctionalGroups(thisMonomerDictionary, thioAlk1Axis, alk1Alk2Axis, thioAACOM, alk1AACOM, alk2AACOM, thioAAIDs, alk1AAIDs, alk2AAIDs)
            atomIDLookupTable[thioID] = ['thio', [x + noAtomsInMolecule + self.noAtomsInMorphology for x in thioAAIDs]]
            atomIDLookupTable[alk1ID] = ['alk1', [x + noAtomsInMolecule + self.noAtomsInMorphology for x in alk1AAIDs]]
            atomIDLookupTable[alk2ID] = ['alk2', [x + noAtomsInMolecule + self.noAtomsInMorphology for x in alk2AAIDs]]
            for thioAAID in thioAAIDs:
                thisMonomerDictionary['unwrapped_position'][thioAAID] = list(np.array(thisMonomerDictionary['unwrapped_position'][thioAAID])+thioTranslation)
                thisMonomerDictionary['velocity'][thioAAID] = self.CGDictionary['velocity'][thioID]
            for alk1AAID in alk1AAIDs:
                thisMonomerDictionary['unwrapped_position'][alk1AAID] = list(np.array(thisMonomerDictionary['unwrapped_position'][alk1AAID])+alk1Translation)
                thisMonomerDictionary['velocity'][alk1AAID] = self.CGDictionary['velocity'][alk1ID]
            for alk2AAID in alk2AAIDs:
                thisMonomerDictionary['unwrapped_position'][alk2AAID] = list(np.array(thisMonomerDictionary['unwrapped_position'][alk2AAID])+alk2Translation)
                thisMonomerDictionary['velocity'][alk2AAID] = self.CGDictionary['velocity'][alk2ID]
            thisMonomerDictionary = helperFunctions.addWrappedPositions(thisMonomerDictionary)
            thisMonomerDictionary = helperFunctions.addMasses(thisMonomerDictionary)
            thisMonomerDictionary = helperFunctions.addDiameters(thisMonomerDictionary)
            
            #### DEBUG CODE ####
            # listOfMasses = []
            # listOfPositions = []
            # for atomNo in thioAAIDs:
            #     listOfPositions.append(thisMonomerDictionary['unwrapped_position'][atomNo])
            #     if ('H' in thisMonomerDictionary['type'][atomNo]):
            #         listOfMasses.append(1.0)
            #     elif ('C' in thisMonomerDictionary['type'][atomNo]):
            #         listOfMasses.append(12.0)
            #     elif ('S' in thisMonomerDictionary['type'][atomNo]):
            #         listOfMasses.append(32.0)
            # print len(listOfPositions), len(listOfMasses)
            # COMPosn = calcCOM(listOfPositions, listOfMasses)
            # print "CGCOM =", thioPosn
            # print "AACOM =", COMPosn
            #####################
                
                                   
            # Now add in the inter-monomer bond.
            # In the template, the inter-monomer bonding carbons are given by type C10 and C1, which are thioAAIDs 0 and 3
            if thioID != moleculeEnds[1]:
                # We are at the start of the molecule, or somewhere in the middle and so an intermonomer bond is required
                # to connect to the next monomer. This bond is between thioAAID 3 and thioAAID 0+25 (again the 25 is hard
                # coded for this template, because it has 25 atoms in it)
                # Create the bond (1 bond)
                thisMonomerDictionary['bond'].append(['C1-C10', 3, 25])
                # Create the angles (4 angles)
                thisMonomerDictionary['angle'].append(['C1-C10-C9', 3, 25, 26])
                thisMonomerDictionary['angle'].append(['C1-C10-S1', 3, 25, 29])
                thisMonomerDictionary['angle'].append(['S1-C1-C10', 4, 3, 25])
                thisMonomerDictionary['angle'].append(['C2-C1-C10', 2, 3, 25])                
                # Create the dihedrals that are necessary to define the forcefield (the commented out ones are not required)
                # thisMonomerDictionary['dihedral'].append(['C1-C10-C9-H1', 3, 25, 26, 49])
                thisMonomerDictionary['dihedral'].append(['C1-C10-C9-C2', 3, 25, 26, 27])
                thisMonomerDictionary['dihedral'].append(['C1-C10-S1-C1', 3, 25, 29, 28])
                # thisMonomerDictionary['dihedral'].append(['C10-C1-S1-C1', 25, 3, 4, 0])
                # thisMonomerDictionary['dihedral'].append(['C10-C1-C2-C9', 25, 3, 2, 1])
                # thisMonomerDictionary['dihedral'].append(['C10-C1-C2-C3', 25, 3, 2, 5])
                # thisMonomerDictionary['dihedral'].append(['S1-C1-C10-C9', 4, 3, 25, 26])
                thisMonomerDictionary['dihedral'].append(['S1-C1-C10-S1', 4, 3, 25, 29])
                # thisMonomerDictionary['dihedral'].append(['C2-C1-C10-C9', 2, 3, 25, 26])
                thisMonomerDictionary['dihedral'].append(['C2-C1-C10-S1', 2, 3, 25, 29])
                
            # Now we update the atom IDs to mirror the fact that we have added an additional monomer to the system
            thisMonomerDictionary = helperFunctions.incrementAtomIDs(thisMonomerDictionary, noAtomsInMolecule)
            noAtomsInMolecule += 25
            AADictionary = self.updateMoleculeDictionary(thisMonomerDictionary, AADictionary)
        AADictionary['natoms'] = noAtomsInMolecule
        for key in ['lx', 'ly', 'lz']:
            AADictionary[key] = thisMonomerDictionary[key]
        AADictionary, terminatingHydrogenIDs = helperFunctions.addTerminatingHydrogens(AADictionary)
        # Update the atomID CG LookupTable with the terminating Hydrogens
        atomIDLookupTable[moleculeEnds[0]][1].append(terminatingHydrogenIDs[0] + self.noAtomsInMorphology)
        atomIDLookupTable[moleculeEnds[1]][1].append(terminatingHydrogenIDs[1] + self.noAtomsInMorphology)
        # Now the molecule is done, we need to add on the correct identifying numbers for all the bonds, angles and dihedrals
        # (just as we did between monomers) for the other molecules in the system, so that they all connect to the right atoms
        AADictionary = helperFunctions.incrementAtomIDs(AADictionary, self.noAtomsInMorphology)
        return AADictionary, atomIDLookupTable


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
    
    def getAATemplatePosition(self, thioAAIDs, alk1AAIDs, alk2AAIDs):
        thioAAPosn = []
        alk1AAPosn = []
        alk2AAPosn = []
        thioAAMasses = []
        alk1AAMasses = []
        alk2AAMasses = []
        # Need to get the masses first, obtained from nist.gov
        for atomID in thioAAIDs:
            if ('C' in self.AATemplateDictionary['type'][atomID]):
                thioAAMasses.append(12.00000)
            elif ('S' in self.AATemplateDictionary['type'][atomID]):
                thioAAMasses.append(31.97207)
            elif ('H' in self.AATemplateDictionary['type'][atomID]):
                thioAAMasses.append(1.00783)
            else:
                raise SystemError('INCORRECT ATOM TYPE', self.AATemplateDictionary['type'][atomID])
            thioAAPosn.append(self.AATemplateDictionary['unwrapped_position'][atomID])
        for atomID in alk1AAIDs:
            if ('C' in self.AATemplateDictionary['type'][atomID]):
                alk1AAMasses.append(12.00000)
            elif ('H' in self.AATemplateDictionary['type'][atomID]):
                alk1AAMasses.append(1.00783)
            else:
                raise SystemError('INCORRECT ATOM TYPE', self.AATemplateDictionary['type'][atomID])
            alk1AAPosn.append(self.AATemplateDictionary['unwrapped_position'][atomID])
        for atomID in alk2AAIDs:
            if ('C' in self.AATemplateDictionary['type'][atomID]):
                alk2AAMasses.append(12.00000)
            elif ('H' in self.AATemplateDictionary['type'][atomID]):
                alk2AAMasses.append(1.00783)
            else:
                raise SystemError('INCORRECT ATOM TYPE', self.AATemplateDictionary['type'][atomID])
            alk2AAPosn.append(self.AATemplateDictionary['unwrapped_position'][atomID])
        # These output as Numpy arrays because we can't do maths with lists
        thiopheneCOM = helperFunctions.calcCOM(thioAAPosn, thioAAMasses)
        alk1COM = helperFunctions.calcCOM(alk1AAPosn, alk1AAMasses)
        alk2COM = helperFunctions.calcCOM(alk2AAPosn, alk2AAMasses)
        # Now move the COM's based on the position of the atoms in self.CGDictionary
        return thiopheneCOM, alk1COM, alk2COM


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
            elif "<velocity" in templateData[lineNo]:
                if (len(self.dataDictionary['velocity']) != 0):
                    for dataToWrite in list(reversed(self.dataDictionary['velocity'])):
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
