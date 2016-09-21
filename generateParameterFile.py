import os
import sys
sys.path.append('./code')
import helperFunctions

def obtainCGMappings(inputDict, templateDict):
    CGSiteTypes = sorted(set(inputDict['type']))
    templateTypes = templateDict['type']
    print "\nThe coarse-grained morphology contains", len(CGSiteTypes), "different coarse-grained sites:"
    for typeName in CGSiteTypes:
        print typeName+" ",
    print "\n\nThe corresponding template file contains", len(templateTypes), "different atoms:"
    for typeName in templateTypes:
        print typeName+" ",
    print "\n"
    CGToTemplateAAIDs = {}
    while True:
        atomQuantityList = []
        for CGSiteType in CGSiteTypes:
            while True:
                try:
                    numberOfAtoms = int(raw_input("How many atoms in the template files are represented by the coarse-grained site, "+str(CGSiteType)+"? "))
                except ValueError:
                    print "Please only enter integers"
                    continue
                if numberOfAtoms <= 0:
                    print "A coarse-grained bead cannot represent zero or fewer atoms"
                else:
                    print ""
                    atomQuantityList.append(numberOfAtoms)
                    break
        if sum(atomQuantityList) > len(templateDict['type']):
            print "The sum of the atoms specified ("+str(sum(atomQuantityList))+") is greater than the number of atoms in the template file ("+str(len(templateDict['type']))+")."
        elif sum(atomQuantityList) < len(templateDict['type']):
            print "The sum of the atoms specified ("+str(sum(atomQuantityList))+") is less than the number of atoms in the template file ("+str(len(templateDict['type']))+")."
        else:
            print ""
            break
    for CGSiteTypeIndex, CGSiteType in enumerate(CGSiteTypes):
        while True:
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
                print ""
                CGToTemplateAAIDs[str(CGSiteType)] = AAIDs
                break
    return CGToTemplateAAIDs


def lookForFiles(inputDirectory, mode):
    while True:
        validFiles = []
        currentDirContents = os.listdir(os.getcwd()+'/'+inputDirectory)
        print "Looking for "+mode+" files..."
        for fileName in currentDirContents:
            if ('.xml' in fileName) and ('template' not in fileName) and ('model' not in fileName):
                validFiles.append(fileName)
        print "\n---=== VALID XML FILES ===---"
        for elementNo in range(len(validFiles)):
            print str(elementNo)+"):", validFiles[elementNo]
        print "\n"+str(elementNo+1)+"): Define alternative path"
        print str(elementNo+2)+"): Exit program\n"
        runThisFile = raw_input("Please pick a file to run (integer, default = 0): ")
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
        useThisDir = raw_input("Please select an output directory: ")
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
        return outputDirectory

def getSigmaVal(morphology):
    while True:
        try:
            sigmaValue = raw_input("Please enter the length scaling value (sigma) for the input morphology, "+str(morphology)+" (float, default = 1.0). Note that this was defined when the CG morphologies were generated and is therefore not necessarily the same as the Lennard-Jones sigma scaling for the atomistic relaxation which will be done automatically: ")
            return float(sigmaValue)
        except:
            print "Please try again."


def checkBonds(inputDict, templateDict, CGToTemplateAAIDs):
    CGToAAIDBonds = {}
    CGBondTypes = sorted(list(set(bond[0] for bond in inputDict['bond'])))
    intraCGBonds = []
    for bond in templateDict['bond']:
        for i, j in enumerate(CGToTemplateAAIDs.values()):
            if bond[1] in j:
                CGSiteAtom1 = i
            if bond[2] in j:
                CGSiteAtom2 = i
        if CGSiteAtom1 != CGSiteAtom2:
            intraCGBonds.append(bond)
    intraCGBonds.sort(key = lambda x: x[1])
    print len(CGBondTypes), "bonds have been detected in the input morphology."
    print len(intraCGBonds), "intra-CG site bonds have been detected in the template."
    print "For each bond in the CG morphology, please indicate the representative atomistic bond in the template:"
    while len(CGBondTypes) > 0:
        print "---=== Defining CG Bond", CGBondTypes[0], "===---"
        while True:
            for index, intraCGBond in enumerate(intraCGBonds):
                print str(index)+"):", intraCGBond
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
                while True:
                    interMonomerBondRaw = raw_input("Please describe the inter-monomer bond in the form 'STR INT(ID of Atom on Monomer 1) INT(ID of Atom on Monomer 2)'. Note that the bond indices will be incremented automatically to reflect which monomer the atom belongs to: ").split(" ")
                    try:
                        # In the code, we want to increment the index of the atom on the second repeat unit so that we know it's a new bond
                        interMonomerBond = [str(interMonomerBondRaw[0]), int(interMonomerBondRaw[1]), int(interMonomerBondRaw[2])+len(templateDict['type'])]
                        break
                    except:
                        print "Please try again."
                CGToAAIDBonds[CGBondTypes[0]] = interMonomerBond
                CGBondTypes.pop(0)
                break
            else:
                CGToAAIDBonds[CGBondTypes[0]] = intraCGBonds[bondSelection]
                CGBondTypes.pop(0)
                intraCGBonds.pop(bondSelection)
                break
    # Add any new constraints required
    if 'interMonomerBond' in locals():
        print "\n The Inter-monomer bond", interMonomerBond, "will probably require additional constraints (angles, dihedrals, impropers)."
        print "Please define any additional constraints now."
        additionalConstraints = []
        while True:
            print "Current Constraints:", additionalConstraints
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
                            while int(interMonomerAngleRaw[i]) > len(templateDict['type']):
                                interMonomerAngleRaw[i] = int(interMonomerAngleRaw[i]) - templateDict['type']
                        interMonomerAngle = [str(interMonomerAngleRaw[0]), int(interMonomerAngleRaw[1]), int(interMonomerAngleRaw[2]), int(interMonomerAngleRaw[3])]
                        print "Angle =", interMonomerAngle
                        # Work out which atoms belong to which monomer for the code
                        for index, atomID in enumerate(interMonomerAngle[1:]):
                            while True:
                                try:
                                    monomerID = raw_input("Does atom with index "+str(atomID)+" in the angle belong to monomer 1 or monomer 2?: ")
                                    interMonomerAngle[index+1] += (int(monomerID)-1)*len(templateDict['type'])
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
                            while int(interMonomerDihedralRaw[i]) > len(templateDict['type']):
                                interMonomerDihedralRaw[i] = int(interMonomerDihedral[i]) - templateDict['type']
                        interMonomerDihedral = [str(interMonomerDihedralRaw[0]), int(interMonomerDihedralRaw[1]), int(interMonomerDihedralRaw[2]), int(interMonomerDihedralRaw[3]), int(interMonomerDihedralRaw[4])]
                        print "Dihedral =", interMonomerDihedral
                        # Work out which atoms belong to which monomer for the code
                        for index, atomID in enumerate(interMonomerDihedral[1:]):
                            while True:
                                try:
                                    monomerID = raw_input("Does atom with index "+str(atomID)+" in the dihedral belong to monomer 1 or monomer 2?: ")
                                    interMonomerDihedral[index+1] += (int(monomerID)-1)*len(templateDict['type'])
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
                            while int(interMonomerImproperRaw[i]) > len(templateDict['type']):
                                interMonomerImproperRaw[i] = int(interMonomerImproper[i]) - templateDict['type']
                        interMonomerImproper = [str(interMonomerImproperRaw[0]), int(interMonomerImproperRaw[1]), int(interMonomerImproperRaw[2]), int(interMonomerImproperRaw[3]), int(interMonomerImproperRaw[4])]
                        print "Improper =", interMonomerImproper
                        # Work out which atoms belong to which monomer for the code
                        for index, atomID in enumerate(interMonomerImproper[1:]):
                            while True:
                                try:
                                    monomerID = raw_input("Does atom with index "+str(atomID)+" in the improper belong to monomer 1 or monomer 2?: ")
                                    interMonomerImproper[index+1] += (int(monomerID)-1)*len(templateDict['type'])
                                    break
                                except:
                                    print "Please try again."
                        additionalConstraints.append(['IMPROPER']+interMonomerImproper)
                        break
                    except:
                        print "Please try again."
            else:
                break
        additionalConstraints.append(interMonomerBond)
        return CGToAAIDBonds, additionalConstraints


def getForcefieldParameters(templateDict, additionalConstraints, manual=True):
    forcefieldParameters = {'BONDCOEFFS':[], 'ANGLECOEFFS':[], 'DIHEDRALCOEFFS':[], 'IMPROPERCOEFFS':[], 'LJPAIRCOEFFS':[], 'DPDPAIRCOEFFS':[]}
    sortedAdditionalConstraints = [[],[],[],[]] # [[]]*4 doesn't work because it just copies the memory address!
    for constraint in additionalConstraints:
        if constraint[0] != 'IMPROPER':
            sortedAdditionalConstraints[len(constraint)-2].append(constraint)
        else:
            sortedAdditionalConstraints[3].append(constraint[1:])
    # Sort out bonds first
    bondTypes = sorted(list(set(bond[0] for bond in templateDict['bond'])) + list(additionalBond[0] for additionalBond in sortedAdditionalConstraints[0]))
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
        angleTypes = sorted(list(set(angle[0] for angle in templateDict['angle'])) + list(additionalAngle[0] for additionalAngle in sortedAdditionalConstraints[1]))
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
    dihedralTypes = sorted(list(set(dihedral[0] for dihedral in templateDict['dihedral'])) + list(additionalDihedral[0] for additionalDihedral in sortedAdditionalConstraints[2]))
    for dihedralType in dihedralTypes:
        if manual == True:
            while True:
                try:
                    dihedralCoeffsRaw = raw_input("Please enter the (unscaled!) V0, V1, V2 and V3 values for the "+str(dihedralType)+" dihedral, separated by a space: ").split(" ")
                    forcefieldParameters['DIHEDRALCOEFFS'].append([dihedralType] + list(map(float, dihedralCoeffsRaw)))
                    break
                except:
                    print "Please try again."
        else:
            forcefieldParameters['DIHEDRALCOEFFS'].append([dihedralType, 1.0, 1.0, 1.0, 1.0])
    # Sort out impropers
    improperTypes = sorted(list(set(improper[0] for improper in templateDict['improper'])) + list(additionalImproper[0] for additionalImproper in sortedAdditionalConstraints[3]))
    for improperType in improperTypes:
        if manual == True:
            while True:
                try:
                    improperCoeffsRaw = raw_input("Please enter the (unscaled!) V0, V1, V2 and V3 values for the "+str(improperType)+" improper, separated by a space: ").split(" ")
                    forcefieldParameters['IMPROPERCOEFFS'].append([improperType] + list(map(float, improperCoeffsRaw)))
                    break
                except:
                    print "Please try again."
        else:
            forcefieldParameters['IMPROPERCOEFFS'].append([improperType, 1.0, 1.0, 1.0, 1.0])
    # Sort out LJ nonbonded pairs (will also add default-valued DPD pair potentials)
    atomTypes = sorted(list(set(templateDict['type'])))
    for atomType1 in atomTypes:
        for atomType2 in atomTypes:
            pairType = str(atomType1)+"-"+str(atomType2)
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


def writeParFile(fileName, INPUTDIR, INPUTMORPHOLOGY, OUTPUTDIR, TEMPLATEDIR, AATEMPLATEFILE, CGToIDDictionary, CGToBondDictionary, forcefieldParameters):
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
        if 'CGTOTEMPLATEAAIDS' in line:
            AAIDLines = ""
            for dictKey in sorted(CGToIDDictionary.keys()):
                AAIDLines += repr(dictKey)+":"+str(CGToIDDictionary[dictKey])+",\\\n"
            parTemplate[lineNo] = AAIDLines
        if 'CGTOTEMPLATEBONDS' in line:
            AABondLines = ""
            for dictKey in sorted(CGToBondDictionary.keys()):
                AABondLines += repr(dictKey)+":"+str(CGToBondDictionary[dictKey])+",\\\n"
            parTemplate[lineNo] = AABondLines
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
    templateDir, AATemplateFile = lookForFiles('templates', 'atomistic template')
    templateDict = helperFunctions.loadMorphologyXML(templateDir+'/'+AATemplateFile)
    # Then, work out the CG mappings
    CGToTemplateAAIDs = obtainCGMappings(inputDict, templateDict)
    # Now, check that the template contains all of the bonds that are specified in the CG morphology
    CGToAAIDBonds, additionalConstraints = checkBonds(inputDict, templateDict, CGToTemplateAAIDs)
    # Next up are the force-field coefficients for all of the bonds, angles and dihedrals specified here
    forcefieldParameters = getForcefieldParameters(templateDict, additionalConstraints, manual=False)
    # Finally, create the parameter file.
    print "The majority of the input parameters are now set."
    fileName = getParFileName()
    print "These, along with some default parameters will be written to "+str(fileName)+"."
    print "You are encouraged to modify "+str(fileName)+" directly now to avoid this process, and also to check that the default variables are desirable."
    writeParFile(fileName, inputDir, inputMorphology, outputDir, templateDir, AATemplateFile, CGToTemplateAAIDs, CGToAAIDBonds, forcefieldParameters)
    print "Parameters written to "+str(fileName)+"."
