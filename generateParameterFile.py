import os
import sys
sys.path.append('./code')
import helperFunctions

def loadMultipleTemplates():
    while True:
        numberOfTemplatesToLoadRaw = raw_input("How many template files would you like to load in order to completely describe the morphology in "+inputDir+"/"+inputMorphology+"? (integer, default = 1): ")
        if len(numberOfTemplatesToLoadRaw) == 0:
            numberOfTemplatesToLoad = 1
            break
        try:
            numberOfTemplatesToLoad = int(numberOfTemplatesToLoadRaw)
            break
        except:
            print "Please try again."
    templateDirs = []
    AATemplateFiles = []
    templateDicts = []
    for i in range(numberOfTemplatesToLoad):
        templateDir, AATemplateFile = lookForFiles('templates', 'atomistic template', ignoreFiles = AATemplateFiles)
        templateDirs.append(templateDir)
        AATemplateFiles.append(AATemplateFile)
        templateDicts.append(helperFunctions.loadMorphologyXML(templateDir+'/'+AATemplateFile))
    return templateDirs, AATemplateFiles, templateDicts


def removeDuplicates(templateDicts):
    allAtomTypes = []
    for templateNo, templateDict in enumerate(templateDicts):
        allAtomTypes.append([])
        allAtomTypes[-1] = list(set(templateDict['type']))
    newAtomTypeList = allAtomTypes[0]
    for templateNo, atomTypesList in enumerate(allAtomTypes):
        if templateNo == 0:
            continue
        for atomType in atomTypesList:
            if atomType in newAtomTypeList:
                # We have a duplicate so fix it
                for i, j in enumerate(atomType):
                    try:
                        int(j)
                        break
                    except:
                        continue
                element = atomType[:i]
                number = int(atomType[i:])
                # Increment the number to find a new atomType that hasn't been used yet
                while True:
                    number += 1
                    newAtomType = element+str(number)
                    if newAtomType in newAtomTypeList:
                        continue
                    # Now we have an AtomType that hasn't been used yet, update the dictionary
                    templateDicts = updateAtomType(templateDicts, templateNo, atomType, newAtomType)
                    newAtomTypeList.append(newAtomType)
                    break
    afterAllAtomTypes = []
    for templateNo, templateDict in enumerate(templateDicts):
        afterAllAtomTypes.append([])
        afterAllAtomTypes[-1] = list(set(templateDict['type']))
    return templateDicts


def updateAtomType(templateDicts, templateNo, oldAtomType, newAtomType):
    atomIDsThatNeedChanging = []
    for atomID, atomType in enumerate(templateDicts[templateNo]['type']):
        if atomType == oldAtomType:
            templateDicts[templateNo]['type'][atomID] = newAtomType
            atomIDsThatNeedChanging.append(atomID)
    for constraintType in ['bond', 'angle', 'dihedral', 'improper']:
        for constraintNo, constraint in enumerate(templateDicts[templateNo][constraintType]):
            for index, atomID in enumerate(constraint):
                if atomID in atomIDsThatNeedChanging:
                    newConstraint0 = constraint[0].split('-')
                    # Sanity check
                    if newConstraint0[index-1] != oldAtomType:
                        raise SystemError('Failed sanity check (see code for details)')
                    newConstraint0[index-1] = newAtomType
                    templateDicts[templateNo][constraintType][constraintNo][0] = '-'.join(newConstraint0)
    return templateDicts


def obtainCGMappings(inputDict, templateDicts, AATemplateFiles, AATemplateDirs):
    CGSiteTypes = sorted(set(inputDict['type']))
    templateTypes = [templateDict['type'] for templateDict in templateDicts]
    print "\nThe coarse-grained morphology contains", len(CGSiteTypes), "different coarse-grained sites:"
    for typeName in CGSiteTypes:
        print typeName+" ",
    totalAtomsInTemplates = sum([len(templateType) for templateType in templateTypes])
    print "\n\nThe corresponding template files contain", totalAtomsInTemplates, "different atoms in total:"
    print "--== Template Files ==--"
    for templateIndex, templateTypesList in enumerate(templateTypes):
        print str(templateIndex)+") "+AATemplateFiles[templateIndex]+" [",
        for typeName in templateTypesList:
            print typeName+" ",
        print "]"
    print "\n"
    CGToTemplateAAIDs = {}
    CGToTemplateDirs = {}
    CGToTemplateFiles = {}
    rigidBodySites = {}
    useAllAtoms = {}
    while True:
        atomQuantityList = []
        for CGSiteType in CGSiteTypes:
            if len(AATemplateFiles) > 1:
                while True:
                    try:
                        templateFileIndex = int(raw_input("Which template file index should be used to represent the coarse-grained site, "+str(CGSiteType)+"? "))
                        break
                    except ValueError:
                        print "Please only enter integers"
                        continue
            else:
                templateFileIndex = 0
            CGToTemplateDirs[CGSiteType] = AATemplateDirs[templateFileIndex]
            CGToTemplateFiles[CGSiteType] = AATemplateFiles[templateFileIndex]
            while True:
                numberOfAtoms = raw_input("How many atoms in the template files are represented by the coarse-grained site, "+str(CGSiteType)+" (leave blank to use all atoms in template file) ? ")
                if len(numberOfAtoms) == 0:
                    numberOfAtoms = len(templateTypes[templateFileIndex])
                    print "Using all", numberOfAtoms, "atoms."
                    useAllAtoms[CGSiteType] = True
                else:
                    try:
                        numberOfAtoms = int(numberOfAtoms)
                        useAllAtoms[CGSiteType] = False
                    except ValueError:
                        print "Please only enter integers"
                        continue
                if numberOfAtoms <= 0:
                    print "A coarse-grained bead cannot represent zero or fewer atoms"
                else:
                    print ""
                    atomQuantityList.append(numberOfAtoms)
                    break
        if sum(atomQuantityList) > totalAtomsInTemplates:
            print "The sum of the atoms specified ("+str(sum(atomQuantityList))+") is greater than the number of atoms across all template files ("+str(totalAtomsInTemplates)+")."
        elif sum(atomQuantityList) < totalAtomsInTemplates:
            print "The sum of the atoms specified ("+str(sum(atomQuantityList))+") is less than the number of atoms across all template files ("+str(totalAtomsInTemplates)+")."
        else:
            print ""
            break
    for CGSiteTypeIndex, CGSiteType in enumerate(CGSiteTypes):
        while True:
            if useAllAtoms[CGSiteType] is True:
                AAIDs = [i for i in range(len(atomQuantityList[CGSiteTypeIndex]))]
            else:
                AAString = str(raw_input("Please type the integer atom IDs (as they appear in the template file) that are described by the coarse-grained bead of type "+str(CGSiteType)+", with each ID separated by a space: "))
                try:
                    AAIDs = map(int, AAString.split(' '))
                except ValueError:
                    print "Syntax error, please try again."
                    continue
            if len(AAIDs) > atomQuantityList[CGSiteTypeIndex]:
                print "The length of the AAIDs specified ("+str(len(AAIDs))+") is greater than the quantity specified previously ("+str(atomQuantityList[CGSiteTypeIndex])+")."
            elif len(AAIDs) < atomQuantityList[CGSiteTypeIndex]:
                print "The length of the AAIDs specified ("+str(len(AAIDs))+") is less than the quantity specified previously ("+str(atomQuantityList[CGSiteTypeIndex])+")."
            else:
                CGToTemplateAAIDs[str(CGSiteType)] = AAIDs
                break
        while True:
            rigidString = str(raw_input("Is site "+str(CGSiteType)+" a rigid body? Leave empty for no or enter the integer atom IDs for the atoms in this group that belong to the rigid body: "))
            if len(rigidString) == 0:
                break
            try:
                rigidIDs = map(int, rigidString.split(' '))
            except ValueError:
                print "Syntax error, please try again."
            for atomID in rigidIDs:
                if atomID not in CGToTemplateAAIDs[str(CGSiteType)]:
                    print "Unexpected atom index found:", atomID
                    continue
            rigidBodySites[str(CGSiteType)] = rigidIDs
            break
        print ""
    return CGToTemplateDirs, CGToTemplateFiles, CGToTemplateAAIDs, rigidBodySites


def lookForFiles(inputDirectory, mode, ignoreFiles=[]):
    while True:
        validFiles = []
        currentDirContents = os.listdir(os.getcwd()+'/'+inputDirectory)
        print "Looking for "+mode+" files..."
        for fileName in currentDirContents:
            if ('.xml' in fileName) and ('template' not in fileName) and ('model' not in fileName) and (fileName not in ignoreFiles):
                validFiles.append(fileName)
        print "\n---=== VALID XML FILES ===---"
        for elementNo in range(len(validFiles)):
            print str(elementNo)+"):", validFiles[elementNo]
        print "\n"+str(elementNo+1)+"): Define alternative path"
        print str(elementNo+2)+"): Exit program\n"
        runThisFile = raw_input("Please pick a file (integer, default = 0): ")
        if len(runThisFile) == 0:
            runThisFile = 0
        else:
            try:
                runThisFile = int(runThisFile)
            except:
                print "Please enter an integer between 0 and", len(validFiles)
                continue
        if (runThisFile < 0) or (runThisFile > len(validFiles)+1):
            print "Please enter an integer between 0 and", len(validFiles)
            continue
        elif runThisFile == len(validFiles):
            inputDirectory = str(raw_input("Please enter the new path where the "+mode+" file is located: "))
            continue
        elif runThisFile == len(validFiles)+1:
            print "Exiting Program..."
            exit()
        return os.getcwd()+'/'+inputDirectory, validFiles[runThisFile]


def lookForDirectory():
    while True:
        potentialDirs = []
        currentDirContents = os.listdir(os.getcwd())
        print "Looking for a suitable output directory..."
        for fileName in currentDirContents:
            if ('output' in fileName) or ('Output' in fileName):
                potentialDirs.append(fileName)
        print "\n---=== VALID OUTPUT DIRECTORIES ===---"
        for elementNo in range(len(potentialDirs)):
            print str(elementNo)+"):", potentialDirs[elementNo]
        print "\n"+str(elementNo+1)+"): Define alternative path"
        print str(elementNo+2)+"): Exit program\n"
        useThisDir = raw_input("Please select an output directory (integer, default = 0): ")
        if len(useThisDir) == 0:
            useThisDir = 0
        else:
            try:
                useThisDir = int(useThisDir)
            except:
                print "Please enter an integer between 0 and", len(potentialDirs)
                continue
        if (useThisDir < 0) or (useThisDir > len(potentialDirs)+1):
            print "Please enter an integer between 0 and", len(potentialDirs)
            continue
        elif useThisDir == len(potentialDirs):
            outputDirectory = str(raw_input("Please enter the new path to the output directory: "))
        elif useThisDir == len(potentialDirs)+1:
            print "Exiting Program..."
            exit()
        else:
            outputDirectory = potentialDirs[useThisDir]
        return os.getcwd()+'/'+outputDirectory

def getSigmaVal(morphology):
    while True:
        try:
            sigmaValue = raw_input("Please enter the length scaling value (sigma) for the input morphology, "+str(morphology)+" (float, default = 1.0). Note that this was defined when the CG morphologies were generated and is therefore not necessarily the same as the Lennard-Jones sigma scaling for the atomistic relaxation which will be done automatically: ")
            return float(sigmaValue)
        except:
            print "Please try again."


def checkBonds(inputDict, templateDicts, CGToTemplateAAIDs, AATemplateFiles):
    CGToAAIDBonds = {}
    CGBondTypes = sorted(list(set(bond[0] for bond in inputDict['bond'])))
    intraCGBonds = []
    intraCGBondTemplates = []
    # First, check if any of the ATOMIDs cross over between the different CGToTemplateAAIDs.values().
    # If not, then CG sites belong to the same template file so we can look for inter-monomer bonds.
    # Otherwise, we have repeating AAIDs, and so we're looking across multiple templates.
    ignoreTheseCGSites = []
    previousAtomIDs = []
    for CGSite in CGToTemplateAAIDs:
        currentSet = set(CGToTemplateAAIDs[CGSite])
        previousSet = set(previousAtomIDs)
        if len(previousSet.intersection(currentSet)) == 0:
            # All new atoms, so don't ignore this site
            previousAtomIDs += CGToTemplateAAIDs[CGSite]
        else:
            # Repeated atomIDs, therefore ignore this site for the inter-monomer calculation
            ignoreTheseCGSites.append(CGSite)
    for templateNo, templateDict in enumerate(templateDicts):
        for bond in templateDict['bond']:
            for i, j in enumerate(CGToTemplateAAIDs.values()):
                if CGToTemplateAAIDs.keys()[i] in ignoreTheseCGSites:
                    continue
                if bond[1] in j:
                    CGSiteAtom1 = i
                if bond[2] in j:
                    CGSiteAtom2 = i
            if CGSiteAtom1 != CGSiteAtom2:
                intraCGBonds.append(bond)
                intraCGBondTemplates.append(AATemplateFiles[templateNo])
    intraCGBonds.sort(key = lambda x: x[2])
    print len(CGBondTypes), "bonds have been detected in the input morphology."
    print len(intraCGBonds), "intra-CG site bonds have been detected in the template."
    while len(CGBondTypes) > 0:
        print "For each bond in the CG morphology, please indicate the representative atomistic bond in the template:"
        print "---=== Defining CG Bond", CGBondTypes[0], "===---"
        while True:
            for index, intraCGBond in enumerate(intraCGBonds):
                print str(index)+"):", intraCGBond, "["+str(intraCGBondTemplates[index])+"]"
            print str(index+1)+"): New Inter-Monomer bond\n"
            bondSelection = raw_input("Please select an option (integer, default = 0): ")
            if len(bondSelection) == 0:
                bondSelection = 0
            else:
                try:
                    bondSelection = int(bondSelection)
                except:
                    print "Please enter an integer between 0 and", len(bondSelection)
                    continue
            if (bondSelection < 0) or (bondSelection > len(intraCGBonds)):
                print "Please enter an integer between 0 and", len(bondSelection)
                continue
            elif bondSelection == len(intraCGBonds):
                if len(AATemplateFiles) > 1:
                    while True:
                        print "--== Template Files ==--"
                        for templateIndex, templateName in enumerate(AATemplateFiles):
                            print str(templateIndex)+") "+templateName
                        try:
                            templateIndexToUse = int(raw_input("Which template file should be used? "))
                            break
                        except ValueError:
                            print "Please only enter integers"
                            continue
                else:
                    templateIndexToUse = 0
                while True:
                    interMonomerBondRaw = raw_input("Please describe the inter-monomer bond in the form 'STR INT(ID of Atom on Monomer 1) INT(ID of Atom on Monomer 2)'. Note that the bond indices will be incremented automatically to reflect which monomer the atom belongs to: ").split(" ")
                    try:
                        # In the code, we want to increment the index of the atom on the second repeat unit so that we know it's a new bond
                        interMonomerBond = [str(interMonomerBondRaw[0]), int(interMonomerBondRaw[1]), int(interMonomerBondRaw[2])+len(templateDicts[templateIndexToUse]['type'])]
                        break
                    except:
                        print "Please try again."
                #CGToAAIDBonds[CGBondTypes[0]] = interMonomerBond - Don't put this here, it will be in `additionalConstraints' instead for ease of implementation in the finegrainer.
                CGBondTypes.pop(0)
                break
            else:
                CGToAAIDBonds[CGBondTypes[0]] = intraCGBonds[bondSelection]
                CGBondTypes.pop(0)
                intraCGBonds.pop(bondSelection)
                break
    # Add any new constraints required
    additionalConstraints = []
    if 'interMonomerBond' in locals():
        print "\n The Inter-monomer bond", interMonomerBond, "will probably require additional constraints (angles, dihedrals, impropers)."
        print "Please define any additional constraints now."
        additionalConstraints.append(interMonomerBond)
        while True:
            print "---=== Defining additional constraints for inter-monomer bond", interMonomerBond, "===---"
            print "0) New Angle Constraint"
            print "1) New Dihedral Constraint"
            print "2) New Improper Constraint"
            print "3) No More New Constraints"
            constraintSelection = raw_input("Please select an option (integer, default = 3): ")
            if len(constraintSelection) == 0:
                constraintSelection = 3
            else:
                try:
                    constraintSelection = int(constraintSelection)
                except:
                    print "Please enter an integer between 0 and 3"
                    continue
            if (constraintSelection < 0) or (constraintSelection > 3):
                print "Please enter an integer between 0 and 3"
                continue
            elif constraintSelection == 0:
                while True:
                    interMonomerAngleRaw = raw_input("Please describe the inter-monomer angle in the form 'STR INT INT INT': ").split(" ")
                    try:
                        # Put a catch in in case the user starts doing the index incrementation manually
                        for i in range(1,4):
                            while int(interMonomerAngleRaw[i]) > len(templateDicts[templateIndexToUse]['type']):
                                interMonomerAngleRaw[i] = int(interMonomerAngleRaw[i]) - templateDicts[templateIndexToUse]['type']
                        interMonomerAngle = [str(interMonomerAngleRaw[0]), int(interMonomerAngleRaw[1]), int(interMonomerAngleRaw[2]), int(interMonomerAngleRaw[3])]
                        print "Angle =", interMonomerAngle
                        # Work out which atoms belong to which monomer for the code
                        for index, atomID in enumerate(interMonomerAngle[1:]):
                            while True:
                                try:
                                    monomerID = raw_input("Does atom with index "+str(atomID)+" in the angle belong to monomer 1 or monomer 2?: ")
                                    interMonomerAngle[index+1] += (int(monomerID)-1)*len(templateDicts[templateIndexToUse]['type'])
                                    break
                                except:
                                    print "Please try again."
                        additionalConstraints.append(interMonomerAngle)
                        break
                    except:
                        print "Please try again."

            elif constraintSelection == 1:
                while True:
                    interMonomerDihedralRaw = raw_input("Please describe the inter-monomer dihedral in the form 'STR INT INT INT INT': ").split(" ")
                    try:
                        # Put a catch in in case the user starts doing the index incrementation manually
                        for i in range(1,5):
                            while int(interMonomerDihedralRaw[i]) > len(templateDicts[templateIndexToUse]['type']):
                                interMonomerDihedralRaw[i] = int(interMonomerDihedral[i]) - templateDicts[templateIndexToUse]['type']
                        interMonomerDihedral = [str(interMonomerDihedralRaw[0]), int(interMonomerDihedralRaw[1]), int(interMonomerDihedralRaw[2]), int(interMonomerDihedralRaw[3]), int(interMonomerDihedralRaw[4])]
                        print "Dihedral =", interMonomerDihedral
                        # Work out which atoms belong to which monomer for the code
                        for index, atomID in enumerate(interMonomerDihedral[1:]):
                            while True:
                                try:
                                    monomerID = raw_input("Does atom with index "+str(atomID)+" in the dihedral belong to monomer 1 or monomer 2?: ")
                                    interMonomerDihedral[index+1] += (int(monomerID)-1)*len(templateDicts[templateIndexToUse]['type'])
                                    break
                                except:
                                    print "Please try again."
                        additionalConstraints.append(interMonomerDihedral)
                        break
                    except:
                        print "Please try again."

            elif constraintSelection == 2:
                while True:
                    interMonomerImproperRaw = raw_input("Please describe the inter-monomer improper in the form 'STR INT INT INT INT': ").split(" ")
                    try:
                        # Put a catch in in case the user starts doing the index incrementation manually
                        for i in range(1,5):
                            while int(interMonomerImproperRaw[i]) > len(templateDicts[templateIndexToUse]['type']):
                                interMonomerImproperRaw[i] = int(interMonomerImproper[i]) - templateDicts[templateIndexToUse]['type']
                        interMonomerImproper = [str(interMonomerImproperRaw[0]), int(interMonomerImproperRaw[1]), int(interMonomerImproperRaw[2]), int(interMonomerImproperRaw[3]), int(interMonomerImproperRaw[4])]
                        print "Improper =", interMonomerImproper
                        # Work out which atoms belong to which monomer for the code
                        for index, atomID in enumerate(interMonomerImproper[1:]):
                            while True:
                                try:
                                    monomerID = raw_input("Does atom with index "+str(atomID)+" in the improper belong to monomer 1 or monomer 2?: ")
                                    interMonomerImproper[index+1] += (int(monomerID)-1)*len(templateDicts[templateIndexToUse]['type'])
                                    break
                                except:
                                    print "Please try again."
                        additionalConstraints.append(['IMPROPER']+interMonomerImproper)
                        break
                    except:
                        print "Please try again."
            else:
                break
    # Now ask about the terminating groups
    print "NOTE: Currently the code only supports terminating molecules with H1 hydrogen atoms. If a termination group is required, please modify the code to incorporate additional termination template files."
    moleculeTerminatingUnits = []
    if len(AATemplateFiles) == 1:
        # Might need a terminating unit if the template is of a polymer
        # CODE THIS IN IF WE EVER NEED IT
#        while True:
#            print "--== Template Files ==--"
#            for templateIndex, templateName in enumerate(AATemplateFiles):
#                print str(templateIndex)+") "+templateName
#            try:
#                templateIndexToUse = int(raw_input("Which of these template files require a terminating unit because they are polymers? "))
#                break
#            except ValueError:
#                print "Please only enter integers"
#                continue
        for templateIndex, templateName in enumerate(AATemplateFiles):
            while True:
                startAtom = str(raw_input("Please indicate the atom ID at the start of the molecule that the terminating hydrogen is attached to in the", templateName, "template file (leave blank for no terminating group, e.g. in small molecules): "))
                if len(startAtom) == 0:
                    break
                try:
                    moleculeTerminatingUnits.append([templateDict[templateIndex]['type'][int(startAtom)]+'-H1', int(startAtom)])
                except:
                    print "Please try again."
                    continue
                endAtom = str(raw_input("Please indicate the atom ID at the end of the molecule that the terminating hydrogen is attached to: "))
                if len(endAtom) == 0:
                    break
                try:
                    moleculeTerminatingUnits.append([templateDict[templateIndex]['type'][int(endAtom)]+'-H1', int(endAtom)])
                    break
                except:
                    print "Please try again."
                    continue
    return CGToAAIDBonds, additionalConstraints, moleculeTerminatingUnits


def getForcefieldParameters(templateDicts, additionalConstraints, moleculeTerminatingUnits, manual=True):
    forcefieldParameters = {'BONDCOEFFS':[], 'ANGLECOEFFS':[], 'DIHEDRALCOEFFS':[], 'IMPROPERCOEFFS':[], 'LJPAIRCOEFFS':[], 'DPDPAIRCOEFFS':[]}
    # Sort out bonds first
    allBonds = []
    allBonds += [templateDict['bond'] for templateDict in templateDicts]
    # Combine the templates into a single list for iterating
    allBonds = [bond for template in allBonds for bond in template]
    if len(allBonds) != 0:
        bondTypes = sorted(list(set(bond[0] for bond in allBonds)))
        for bondType in bondTypes:
            if manual == True:
                while True:
                    try:
                        bondCoeffsRaw = raw_input("Please enter the (unscaled!) k and r0 values for the "+str(bondType)+" bond, separated by a space: ").split(" ")
                        forcefieldParameters['BONDCOEFFS'].append([bondType] + list(map(float, bondCoeffsRaw)))
                        break
                    except:
                        print "Please try again."
            else:
                forcefieldParameters['BONDCOEFFS'].append([bondType, 1.0, 1.0])
    # Sort out angles
    allAngles = []
    allAngles += [templateDict['angle'] for templateDict in templateDicts]
    # Combine the templates into a single list for iterating
    allAngles = [angle for template in allAngles for angle in template]
    if len(allAngles) != 0:
        angleTypes = sorted(list(set(angle[0] for angle in allAngles)))
        for angleType in angleTypes:
            if manual == True:
                while True:
                    try:
                        angleCoeffsRaw = raw_input("Please enter the (unscaled!) k and theta0 values for the "+str(angleType)+" angle, separated by a space: ").split(" ")
                        forcefieldParameters['ANGLECOEFFS'].append([angleType] + list(map(float, angleCoeffsRaw)))
                        break
                    except:
                        print "Please try again."
            else:
                forcefieldParameters['ANGLECOEFFS'].append([angleType, 1.0, 1.0])
    # Sort out dihedrals
    allDihedrals = []
    allDihedrals += [templateDict['dihedral'] for templateDict in templateDicts]
    # Combine the templates into a single list for iterating
    allDihedrals = [dihedral for template in allDihedrals for dihedral in template]
    if len(allDihedrals) != 0:
        dihedralTypes = sorted(list(set(dihedral[0] for dihedral in allDihedrals)))
        for dihedralType in dihedralTypes:
            if manual == True:
                while True:
                    try:
                        dihedralCoeffsRaw = raw_input("Please enter the (unscaled!) V0, V1, V2, V3 and V4 values for the "+str(dihedralType)+" dihedral, separated by a space: ").split(" ")
                        forcefieldParameters['DIHEDRALCOEFFS'].append([dihedralType] + list(map(float, dihedralCoeffsRaw)))
                        break
                    except:
                        print "Please try again."
            else:
                forcefieldParameters['DIHEDRALCOEFFS'].append([dihedralType, 1.0, 1.0, 1.0, 1.0, 1.0])
    # Sort out impropers
    allImpropers = []
    allImpropers += [templateDict['improper'] for templateDict in templateDicts]
    # Combine the templates into a single list for iterating
    allImpropers = [improper for template in allImpropers for improper in template]
    if len(allImpropers) != 0:
        improperTypes = sorted(list(set(improper[0] for improper in allImpropers)))
        for improperType in improperTypes:
            if manual == True:
                while True:
                    try:
                        improperCoeffsRaw = raw_input("Please enter the (unscaled!) k and chi values for the "+str(improperType)+" improper, separated by a space: ").split(" ")
                        forcefieldParameters['IMPROPERCOEFFS'].append([improperType] + list(map(float, improperCoeffsRaw)))
                        break
                    except:
                        print "Please try again."
            else:
                forcefieldParameters['IMPROPERCOEFFS'].append([improperType, 1.0, 1.0])
    # Sort out LJ nonbonded pairs (will also add default-valued DPD pair potentials)
    allAtomTypes = []
    allAtomTypes += [templateDict['type'] for templateDict in templateDicts]
    # Combine the templates into a single list for iterating
    allAtomTypes = [atomType for template in allAtomTypes for atomType in template]
    atomTypes = sorted(list(set(allAtomTypes)), key=lambda x: helperFunctions.convertStringToInt(x))
    pairTypesIncluded = []
    for atomType1 in atomTypes:
        for atomType2 in atomTypes:
            pairType = str(atomType1)+"-"+str(atomType2)
            pairTypeReversed = str(atomType2)+"-"+str(atomType1)
            if (pairType in pairTypesIncluded) or (pairTypeReversed in pairTypesIncluded):
                continue
            if manual == True:
                while True:
                    try:
                        pairCoeffsRaw = raw_input("Please enter the (unscaled!) epsilon and sigma values for the "+str(pairType)+" pair potential, separated by a space: ").split(" ")
                        forcefieldParameters['LJPAIRCOEFFS'].append([pairType] + list(map(float, pairCoeffsRaw)))
                        forcefieldParameters['DPDPAIRCOEFFS'].append([pairType] + [1.0, 2.96])
                        break
                    except:
                        print "Please try again."
            else:
                forcefieldParameters['LJPAIRCOEFFS'].append([pairType, 1.0, 1.0])
                forcefieldParameters['DPDPAIRCOEFFS'].append([pairType, 1.0, 1.0])
            pairTypesIncluded.append(pairType)
    # Sort out additional constraints
    for constraint in (additionalConstraints + moleculeTerminatingUnits):
        if len(constraint) == 3:
            # Bond
            forcefieldParameters['BONDCOEFFS'].append([constraint[0], 1.0, 1.0])
        elif len(constraint) == 4:
            # Angle
            forcefieldParameters['ANGLECOEFFS'].append([constraint[0], 1.0, 1.0])
        elif len(constraint) == 5:
            # Dihedral
            forcefieldParameters['DIHEDRALCOEFFS'].append([constraint[0], 1.0, 1.0, 1.0, 1.0])
        elif len(constraint) == 6:
            # Improper
            forcefieldParameters['IMPROPERCOEFFS'].append([constraint[1], 1.0, 1.0])
    return forcefieldParameters


def getParFileName():
    filesList = os.listdir(os.getcwd())
    parCounter = 0
    for fileName in filesList:
        if ("par" in fileName) and (".py" in fileName):
            parCounter += 1
    parCounter = str(parCounter)
    while len(parCounter) < 2:
        parCounter = '0'+parCounter
    testFileName = "par"+parCounter+".py"
    while testFileName in filesList:
        testFileName = testFileName[:-4]+"_2.py"
    return testFileName


def writeParFile(fileName, INPUTDIR, INPUTMORPHOLOGY, INPUTSIGMA, OUTPUTDIR, CGToTemplateDirsDictionary, CGToTemplateFilesDictionary, CGToIDDictionary, CGToBondDictionary, forcefieldParameters, additionalInterMonomerConstraints, AAIDsInRigidBodies, moleculeTerminatingUnitConnections):
    with open(os.getcwd()+'/templates/parTemplate.py', 'r') as parTemplateFile:
        parTemplate = parTemplateFile.readlines()
    # Define these in the locals() before we iterate to prevent a dictionary chaning size error
    variableName = None
    variableValue = None
    for lineNo, line in enumerate(parTemplate):
        # Change the simple variables
        for variableName, variableValue in locals().iteritems():
            if variableName in line:
                parTemplate[lineNo] = line.replace(variableName, repr(variableValue))
        # Change the dictionaries
        if 'CGTOTEMPLATEDIRS' in line:
            templateLines = ""
            for dictKey in sorted(CGToTemplateDirsDictionary.keys()):
                templateLines += repr(dictKey)+":"+repr(CGToTemplateDirsDictionary[dictKey])+",\\\n"
            parTemplate[lineNo] = templateLines
        if 'CGTOTEMPLATEFILES' in line:
            templateLines = ""
            for dictKey in sorted(CGToTemplateFilesDictionary.keys()):
                templateLines += repr(dictKey)+":"+repr(CGToTemplateFilesDictionary[dictKey])+",\\\n"
            parTemplate[lineNo] = templateLines
        if 'CGTOTEMPLATEAAIDS' in line:
            AAIDLines = ""
            for dictKey in sorted(CGToIDDictionary.keys()):
                AAIDLines += repr(dictKey)+":"+str(CGToIDDictionary[dictKey])+",\\\n"
            parTemplate[lineNo] = AAIDLines
        elif 'CGTOTEMPLATEBONDS' in line:
            AABondLines = ""
            for dictKey in sorted(CGToBondDictionary.keys()):
                AABondLines += repr(dictKey)+":"+str(CGToBondDictionary[dictKey])+",\\\n"
            parTemplate[lineNo] = AABondLines
        elif 'ADDITIONALCONSTRAINTS' in line:
            addConstLines = ""
            for constraint in additionalInterMonomerConstraints:
                addConstLines += str(constraint)+",\\\n"
            parTemplate[lineNo] = addConstLines
        elif 'RIGIDBODYSITES' in line:
            rigidBodyLines = ""
            for dictKey in sorted(AAIDsInRigidBodies.keys()):
                rigidBodyLines += repr(dictKey)+":"+str(AAIDsInRigidBodies[dictKey])+",\\\n"
            parTemplate[lineNo] = rigidBodyLines
        elif 'TERMINATINGCONNECTIONS' in line:
            terminationLines = ""
            for constraint in moleculeTerminatingUnitConnections:
                terminationLines += str(constraint)+",\\\n"
            parTemplate[lineNo] = terminationLines
        # Change the forcefield lists
        for FFCoeffsName, FFCoeffsValue in forcefieldParameters.iteritems():
            if FFCoeffsName in line:
                constraintLines = ""
                for constraint in FFCoeffsValue:
                    constraintLines += str(constraint)+",\\\n"
                parTemplate[lineNo] = constraintLines
    # Write the file
    with open(os.getcwd()+'/'+fileName, 'w+') as parFile:
        parFile.writelines(parTemplate)


if __name__ == "__main__":
    # First, get the right morphology and template file by asking the user
    inputDir, inputMorphology = lookForFiles('inputCGMorphs', 'input CG morphology')
    inputSigma = getSigmaVal(inputMorphology)
    inputDict = helperFunctions.loadMorphologyXML(inputDir+'/'+inputMorphology, sigma=inputSigma)
    outputDir = lookForDirectory()
    templateDirs, AATemplateFiles, templateDicts = loadMultipleTemplates()
    # If we've loaded two (or more) different templates, we need to make sure there are no overlaps between atom types in the two. Otherwise we will get the wrong bond coefficients and forces etc.
    if len(AATemplateFiles) > 1:
        templateDicts = removeDuplicates(templateDicts)
    # Then, work out the CG mappings
    CGToTemplateDirs, CGToTemplateFiles, CGToTemplateAAIDs, rigidBodySites = obtainCGMappings(inputDict, templateDicts, AATemplateFiles, templateDirs)
    # Now, check that the template contains all of the bonds that are specified in the CG morphology
    CGToAAIDBonds, additionalConstraints, moleculeTerminatingUnits = checkBonds(inputDict, templateDicts, CGToTemplateAAIDs, AATemplateFiles)
    # Next up are the force-field coefficients for all of the bonds, angles and dihedrals specified here
    forcefieldParameters = getForcefieldParameters(templateDicts, additionalConstraints, moleculeTerminatingUnits, manual=False)
    # Finally, create the parameter file.
    print "The majority of the input parameters are now set."
    fileName = getParFileName()
    print "These, along with some default parameters will be written to "+str(fileName)+"."
    print "You are encouraged to modify "+str(fileName)+" directly now to avoid this process, and also to check that the default variables are desirable."
    writeParFile(fileName, inputDir, inputMorphology, inputSigma, outputDir, CGToTemplateDirs, CGToTemplateFiles, CGToTemplateAAIDs, CGToAAIDBonds, forcefieldParameters, additionalConstraints, rigidBodySites, moleculeTerminatingUnits)
    print "Parameters written to "+str(fileName)+"."
