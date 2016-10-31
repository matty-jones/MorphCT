from hoomd_script import *
import numpy as np
import copy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cPickle as pickle
import helperFunctions
import sys

DEBUGWriteDCDFiles = True

class ExitHoomd(Exception):
    def __init__(self, string):
        self.string = string + " At Timestep = " + str(get_step())
    def __str__(self):
        return self.string


class MDPhase:
    def __init__(self, AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, phaseNumber, inputFile, outputFile, sScale, eScale):
        self.AAMorphologyDict = AAMorphologyDict
        self.CGMorphologyDict = CGMorphologyDict
        self.CGToAAIDMaster = CGToAAIDMaster
        self.inputFile = inputFile
        self.outputFile = outputFile
        self.sScale = sScale
        self.eScale = eScale
        print "Note: self.sScale =", self.sScale, "self.eScale =", self.eScale
        self.phaseNumber = phaseNumber
        for key in parameterDict.keys():
            self.__dict__[key] = parameterDict[key]
        for key in ['temperatures', 'taus', 'pairTypes', 'bondTypes', 'angleTypes', 'dihedralTypes', 'integrationTargets', 'timesteps', 'durations', 'terminationConditions', 'groupAnchorings']:
            if self.phaseNumber + 1 > len(parameterDict[key]):
                self.__dict__[key[:-1]] = parameterDict[key][0]
            else:
                self.__dict__[key[:-1]] = parameterDict[key][phaseNumber]
        self.system = init.read_xml(filename = inputFile)
        self.rigidGroup, self.nonRigidGroup = self.getIntegrationGroups()
        self.outputLogFileName = self.outputDir+'/'+self.morphology[:-4]+'/morphology/energies_'+self.morphology[:-4]+'.log'
        self.logQuantities = ['temperature', 'pressure', 'volume', 'potential_energy', 'kinetic_energy', 'bond_'+self.bondType+'_energy', 'angle_'+self.angleType+'_energy', 'dihedral_'+self.dihedralType+'_energy']
        self.getFFCoeffs()

    def optimiseStructure(self):
        if DEBUGWriteDCDFiles is True:
            self.dumpDCD = dump.dcd(filename = self.outputFile.replace('xml', 'dcd'), period = self.duration / 100.0, overwrite = True)
        else:
            self.dumpDCD = None
        self.step = integrate.mode_standard(dt = self.timestep)
        self.rigidInt = integrate.nvt_rigid(group = self.rigidGroup, T = self.temperature, tau = self.tau)
        self.nonRigidInt = integrate.nvt(group = self.nonRigidGroup, T = self.temperature, tau = self.tau)
        if self.phaseNumber == 0:
            logOverwrite = True
        else:
            logOverwrite = False
        self.energyLog = analyze.log(filename = self.outputLogFileName, quantities = self.logQuantities, period = self.duration / 1000.0, overwrite = logOverwrite)
        callback = None
        if self.terminationCondition == 'KEmin':
            self.loadFromSnapshot = False
            self.lowestKE = 9e999
            self.KEIncreased = 0
            self.firstKEValue = True
            callback = analyze.callback(callback = self.checkKE, period = self.duration / 1000.0)
        print "---=== BEGINNING MOLECULAR DYNAMICS PHASE", self.phaseNumber+1, "===---"
        try:
            run(self.duration)
        except ExitHoomd as exitMessage:
            print exitMessage
        if self.terminationCondition == 'KEmin':
            if self.loadFromSnapshot is True:
                print "Loading from snapshot..."
                self.system.restore_snapshot(self.snapshotToLoad)
                del self.snapshotToLoad
        self.dumpXML = dump.xml(filename = self.outputFile, position = True, image = True, type = True, mass = True, diameter = True, body = True, charge = True, bond = True, angle = True, dihedral = True, improper = True)
        self.rigidInt.disable()
        self.nonRigidInt.disable()
        del self.system, self.dumpDCD, self.step, self.rigidInt, self.nonRigidInt, self.rigidGroup, self.nonRigidGroup, callback, self.dumpXML, self.pairClass, self.bondClass, self.angleClass, self.dihedralClass, self.improperClass, self.energyLog
        init.reset()

    def checkKE(self, timestepNumber):
        currentKE = self.energyLog.query('kinetic_energy')
        if self.firstKEValue == False:
            if currentKE >= self.lowestKE:
                if self.KEIncreased == 5:
                    # Found the lowest KE point for at least 5 timesteps
                    del self.firstKEValue, self.lowestKE, self.KEIncreased
                    raise ExitHoomd("Lowest Energy Condition Met")
                self.KEIncreased += 1
            else:
                # Maybe at the lowest KE point so store snapshot
                self.KEIncreased = 0
                self.loadFromSnapshot = True
                self.snapshotToLoad = self.system.take_snapshot(all=True)
                self.lowestKE = currentKE
        else:
            self.firstKEValue = False
        return 0

    def getFFCoeffs(self):
        # Set Pair Coeffs
        self.pairClass = None
        if self.pairType != 'none':
            self.logQuantities.append('pair_'+self.pairType+'_energy')
            # Hoomd crashes if you don't specify all pair combinations, so need to make sure we do this.
            atomTypes = sorted(list(set(self.AAMorphologyDict['type'])), key = lambda x: helperFunctions.convertStringToInt(x))
            allPairTypes = []
            for atomType1 in atomTypes:
                for atomType2 in atomTypes:
                    pairType = str(atomType1)+'-'+str(atomType2)
                    reversePairType = str(atomType2)+'-'+str(atomType1)
                    if (pairType not in allPairTypes) and (reversePairType not in allPairTypes):
                        allPairTypes.append(pairType)
            if self.pairType == 'dpd':
                self.pairClass = pair.dpd(r_cut = self.pairRCut*self.sScale, T = self.temperature)
                for pairCoeff in self.dpdCoeffs:
                    self.pairClass.pair_coeff.set(pairCoeff[0].split('-')[0], pairCoeff[0].split('-')[1], A = pairCoeff[1] * self.eScale, r_cut = pairCoeff[2] * self.sScale * 2**(1./6.), gamma = self.pairDPDGammaVal)
                    try:
                        allPairTypes.remove(pairCoeff[0])
                    except:
                        allPairTypes.remove('-'.join(pairCoeff[0].split('-')[::-1]))  # The reverse way round
                # Now allPairTypes is just the unspecified pair potentials, so set these interactions to zero
                for pairType in allPairTypes:
                    self.pairClass.pair_coeff.set(pairType.split('-')[0], pairType.split('-')[1], A = 0.0, r_cut = 0.0, gamma = 0.0)
            elif self.pairType == 'lj':
                self.pairClass = pair.lj(r_cut = self.pairRCut * self.sScale)
                self.pairClass.set_params(mode = 'xplor')
                for pairCoeff in self.ljCoeffs:
                    self.pairClass.pair_coeff.set(pairCoeff[0].split('-')[0], pairCoeff[0].split('-')[1], epsilon = pairCoeff[1] * self.eScale, sigma = pairCoeff[2] * self.sScale)
                    try:
                        allPairTypes.remove(pairCoeff[0])
                    except:
                        allPairTypes.remove('-'.join(pairCoeff[0].split('-')[::-1]))
                # Now allPairTypes is just the unspecified pair potentials, so set these interactions to zero
                for pairType in allPairTypes:
                    self.pairClass.pair_coeff.set(pairType.split('-')[0], pairType.split('-')[1], epsilon = 0.0, sigma = 0.0)
            else:
                raise SystemError('Non-DPD/LJ pair potentials not yet hard-coded! Please describe how to interpret them on this line.')
        # Set Bond Coeffs
        # Real bonds
        if self.bondType == 'harmonic':
            self.bondClass = bond.harmonic()
            for bondCoeff in self.bondCoeffs:
                # [k] = kcal mol^{-1} \AA^{-2} * episilon/sigma^{2}, [r0] = \AA * sigma^{2}
                self.bondClass.bond_coeff.set(bondCoeff[0], k = bondCoeff[1] * (self.eScale / (self.sScale**2)), r0 = bondCoeff[2] * self.sScale)
        # Ghost bonds
        # If there is no anchoring, rather than change the XML, just set the bond k values to 0.
            if self.groupAnchoring == 'all':
                groupAnchoringTypes = ['X'+CGType for CGType in self.CGToTemplateAAIDs.keys()]
            elif self.groupAnchoring == 'none':
                groupAnchoringTypes = []
            else:
                groupAnchoringTypes = ['X'+CGType for CGType in self.groupAnchoring.split(',')]
            anchorBondTypes = []
            noAnchorBondTypes = []
            for bondType in self.AAMorphologyDict['bond']:
                if 'X' in bondType[0]:
                    atomType1 = bondType[0].split('-')[0]
                    atomType2 = bondType[0].split('-')[1]
                    if (atomType1 in groupAnchoringTypes) or (atomType2 in groupAnchoringTypes):
                        if bondType[0] not in anchorBondTypes:
                            anchorBondTypes.append(bondType[0])
                    else:
                        if bondType[0] not in noAnchorBondTypes:
                            noAnchorBondTypes.append(bondType[0])
            for bondType in anchorBondTypes:
                self.bondClass.bond_coeff.set(bondType, k = 1E6, r0 = 0)
            for bondType in noAnchorBondTypes:
                self.bondClass.bond_coeff.set(bondType, k = 0, r0 = 0)
        else:
            raise SystemError('Non-harmonic bond potentials not yet hard-coded! Please describe how to interpret them on this line.')
        # Set Angle Coeffs
        if self.angleType == 'harmonic':
            self.angleClass = angle.harmonic()
            for angleCoeff in self.angleCoeffs:
                # [k] = kcal mol^{-1} rad^{-2} * epsilon, [t] = rad
                self.angleClass.set_coeff(angleCoeff[0], k = angleCoeff[1] * self.eScale, t0 = angleCoeff[2])
        else:
            raise SystemError('Non-harmonic angle potentials not yet hard-coded! Please describe how to interpret them on this line.')
        # Set Dihedral Coeffs
        if self.dihedralType == 'table':
            self.dihedralClass = dihedral.table(width = 1000)
            for dihedralCoeff in self.dihedralCoeffs:
                self.dihedralClass.dihedral_coeff.set(dihedralCoeff[0], func = multiHarmonicTorsion, coeff = dict(V0 = dihedralCoeff[1] * self.eScale, V1 = dihedralCoeff[2] * self.eScale, V2 = dihedralCoeff[3] * self.eScale, V3 = dihedralCoeff[4] * self.eScale, V4 = dihedralCoeff[5] * self.eScale))
        elif self.dihedralType == 'opls':
            self.dihedralClass = dihedral.opls()
            for dihedralCoeff in self.dihedralCoeffs:
                self.dihedralClass.set_coeff(dihedralCoeff[0], k1 = dihedralCoeff[1] * self.eScale, k2 = dihedralCoeff[2] * self.eScale, k3 = dihedralCoeff[3] * self.eScale, k4 = dihedralCoeff[4] * self.eScale)
        else:
            raise SystemError('Non-tabulated dihedral potentials not yet hard-coded! Please describe how to interpret them on this line.')
        # Set Improper Coeffs
        self.improperClass = None
        if len(self.improperCoeffs) > 0:
            self.improperClass = improper.harmonic()
            for improperCoeff in self.improperCoeffs:
                self.improperClass.improper_coeff.set(improperCoeff[0], k = improperCoeff[1] * self.eScale, chi = improperCoeff[2])

    def getIntegrationGroups(self):
        # Based on input parameter, return all non-rigid and rigid atoms to be integrated over
        if self.integrationTarget == 'all':
            integrationTypes = self.CGToTemplateAAIDs.keys()
        else:
            integrationTypes = self.integrationTarget.split(',')
        # Add in any rigid ghost particles that might need to be integrated too
        ghostIntegrationTypes = ['R' + typeName for typeName in integrationTypes]
        atomIDsToIntegrate = []
        for molecule in self.CGToAAIDMaster:
            for CGSiteID in molecule.keys():
                if molecule[CGSiteID][0] in integrationTypes:
                    atomIDsToIntegrate += molecule[CGSiteID][1]
        for atomID, atomType in enumerate(self.AAMorphologyDict['type']):
            if atomType in ghostIntegrationTypes:
                atomIDsToIntegrate.append(atomID)
        integrateGroup = group.tag_list(name = "integrateGroup", tags = atomIDsToIntegrate)
        rigidGroup = group.intersection(name = "rigidGroup", a = group.rigid(), b = integrateGroup)
        nonRigidGroup = group.difference(name = "nonRigidGroup", a = integrateGroup, b = rigidGroup)
        return rigidGroup, nonRigidGroup


def multiHarmonicTorsion(theta, V0, V1, V2, V3, V4):
    V = V0 + V1*np.cos(theta) + V2*((np.cos(theta))**2) + V3*((np.cos(theta))**3) + V4*((np.cos(theta))**4)
    F = V1*np.sin(theta) + 2*V2*np.cos(theta)*np.sin(theta) + 3*V3*((np.cos(theta))**2)*np.sin(theta) + 4*V4*((np.cos(theta))**3)*np.sin(theta)
    return (V, F)


def obtainScaleFactors(parameterDict):
    print "Obtaining correct scaling for epsilon and sigma..."
    # First determine the correct length scaling
    largestSigma = max(map(float, np.array(parameterDict['ljCoeffs'])[:,2]))
    largestEpsilon = max(map(float, np.array(parameterDict['ljCoeffs'])[:,1]))
    initialMorphology = helperFunctions.loadMorphologyXML(parameterDict['outputDir']+'/'+parameterDict['morphology'][:-4]+'/morphology/'+parameterDict['morphology'])
    return 1/float(largestSigma), 1/float(largestEpsilon)


def scaleMorphology(initialMorphology, parameterDict, sScale, eScale):
    print "Scaling morphology by sigma =", str(1/sScale)+"..."
    if sScale != 1.0:
        helperFunctions.scale(initialMorphology, sScale)
    helperFunctions.writeMorphologyXML(initialMorphology, parameterDict['outputDir']+'/'+parameterDict['morphology'][:-4]+'/morphology/phase0_'+parameterDict['morphology'])


def execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict):
    currentFiles = os.listdir(parameterDict['outputDir']+'/'+parameterDict['morphology'][:-4]+'/morphology')
    sScale, eScale = obtainScaleFactors(parameterDict)
    if (parameterDict['overwriteCurrentData'] is False) and ('phase0_'+parameterDict['morphology'] in currentFiles):
        pass
    else:
        scaleMorphology(AAMorphologyDict, parameterDict, sScale, eScale)
    # Reset logfile
    try:
        os.remove(parameterDict['outputDir']+'/'+parameterDict['morphology'][:-4]+'/morphology/energies_'+parameterDict['morphology'][:-4]+'.log')
    except OSError:
        pass
    for phaseNo in range(parameterDict['numberOfPhases']):
        inputFile = 'phase'+str(phaseNo)+'_'+parameterDict['morphology']
        outputFile = 'phase'+str(phaseNo+1)+'_'+parameterDict['morphology']
        if outputFile in currentFiles:
            if parameterDict['overwriteCurrentData'] is False:
                print outputFile, "already exists. Skipping..."
                continue
        MDPhase(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, phaseNo, parameterDict['outputDir']+'/'+parameterDict['morphology'][:-4]+'/morphology/'+inputFile, parameterDict['outputDir']+'/'+parameterDict['morphology'][:-4]+'/morphology/'+outputFile, sScale, eScale).optimiseStructure()
    # Now all phases are complete, we need to remove the ghost particles from the system
    print "Removing ghost particles to create final output..."
    removeGhostParticles(parameterDict['outputDir']+'/'+parameterDict['morphology'][:-4]+'/morphology/'+outputFile, parameterDict['outputDir']+'/'+parameterDict['morphology'][:-4]+'/morphology/final_'+parameterDict['morphology'])
    return AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict


def removeGhostParticles(lastPhaseXML, outputFileName):
    finalMorphology = helperFunctions.loadMorphologyXML(lastPhaseXML)
    atomIDsToRemove = []
    for atomID, atomType in enumerate(finalMorphology['type']):
        if (atomType[0] == 'X') or (atomType[0] == 'R'):
            # This is a ghost particle
            atomIDsToRemove.append(atomID)
    atomIDsToRemove.sort(reverse = True)
    atomAttribs = ['position', 'image', 'type', 'mass', 'diameter', 'body', 'charge']
    for atomID in atomIDsToRemove:
        for key in atomAttribs:
            finalMorphology[key].pop(atomID)
    finalMorphology['natoms'] -= len(atomIDsToRemove)
    atomConstraints = ['bond', 'angle', 'dihedral', 'improper']
    for key in atomConstraints:
        constraintsToRemove = []
        for constraintNo, constraint in enumerate(finalMorphology[key]):
            for atomID in constraint[1:]:
                if (atomID in atomIDsToRemove) and (constraintNo not in constraintsToRemove):
                    constraintsToRemove.append(constraintNo)
        constraintsToRemove.sort(reverse = True)
        for constraintNo in constraintsToRemove:
            finalMorphology[key].pop(constraintNo)
    helperFunctions.writeMorphologyXML(finalMorphology, outputFileName)


def checkSaveDirectory(directory):
    saveDirectoryFiles = os.listdir(directory+'/morphology')
    runPhase1 = True
    runPhase2 = True
    runPhase3 = True
    runPhase4 = True
    runPhase5 = True
    runPhase6 = True
    runPhase7 = True
    continuePhase7 = False
    continueFile = None
    for fileName in saveDirectoryFiles:
        if ('relaxed' in fileName) and ('xml' in fileName):
            print "Calculations already complete for this morphology."
            return False, False, False, False, False
        elif ('phase1' in fileName) and ('xml' in fileName):
            runPhase1 = False
        elif ('phase2' in fileName) and ('xml' in fileName):
            runPhase2 = False
        elif ('phase3' in fileName) and ('xml' in fileName):
            runPhase3 = False
        elif ('phase4' in fileName) and ('xml' in fileName):
            runPhase4 = False
        elif ('phase5' in fileName) and ('xml' in fileName):
            runPhase5 = False
        elif ('phase6' in fileName) and ('xml' in fileName):
            runPhase6 = False
        elif ('temp' in fileName) and ('xml' in fileName):
            runPhase7 = False
            continuePhase7 = True
            continueFile = saveDirectory+'/'+fileName
    return [runPhase1, runPhase2, runPhase3, runPhase4, runPhase5, runPhase6, runPhase7, continuePhase7, continueFile]

if __name__ == "__main__":
    try:
        pickleFile = sys.argv[1]
    except:
        print "Please specify the pickle file to load to continue the pipeline from this point."
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict = helperFunctions.loadPickle(pickleFile)
    execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict)
