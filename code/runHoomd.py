from hoomd_script import *
import numpy as np
import helperFunctions
import sys

class ExitHoomd(Exception):
    '''This class is raised to terminate a HOOMD simulation mid-run for a particular reason (e.g. minimum KE found)'''
    def __init__(self, string):
        self.string = string + " At Timestep = " + str(get_step())

    def __str__(self):
        return self.string


class MDPhase:
    def __init__(self, AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, phaseNumber, inputFile, outputFile, sScale, eScale):
        '''The MDPhase class respresents a single MD simulation using the given parameters'''
        self.AAMorphologyDict = AAMorphologyDict
        self.CGMorphologyDict = CGMorphologyDict
        self.CGToAAIDMaster = CGToAAIDMaster
        self.inputFile = inputFile
        self.outputFile = outputFile
        self.sScale = sScale
        self.eScale = eScale
        self.phaseNumber = phaseNumber
        # Obtain the parXX.py parameters
        for key in list(parameterDict.keys()):
            self.__dict__[key] = parameterDict[key]
        # Get the phase-specific simulation parameters
        for key in ['temperatures', 'taus', 'pairTypes', 'bondTypes', 'angleTypes', 'dihedralTypes', 'integrationTargets', 'timesteps', 'durations', 'terminationConditions', 'groupAnchorings', 'DCDFileDumpsteps']:
            # If the phase-specific parameter is not defined for this phase number then use the first one.
            if self.phaseNumber + 1 > len(parameterDict[key]):
                self.__dict__[key[:-1]] = parameterDict[key][0]
            # If the phase-specific parameter is specified for this phase then use this parameter
            else:
                self.__dict__[key[:-1]] = parameterDict[key][phaseNumber]
        self.system = self.loadSystem(inputFile)
        # Determine the required groups so we can use the correct integrator each time
        self.rigidGroup, self.nonRigidGroup, self.integrationTypes = self.getIntegrationGroups()
        self.outputLogFileName = self.outputMorphDir + '/' + self.morphology[:-4] + '/morphology/energies_' + self.morphology[:-4] + '.log'
        # Determine which quantities should be logged during the simulation phase
        self.logQuantities = ['temperature', 'pressure', 'volume', 'potential_energy', 'kinetic_energy', 'bond_' + self.bondType + '_energy', 'angle_' + self.angleType + '_energy', 'dihedral_' + self.dihedralType + '_energy']
        # Set the bond coefficients
        self.getFFCoeffs()

    def loadSystem(self, inputFile):
        context.initialize()
        # Load the previous phases' xml for continuation
        systemXML = init.read_xml(filename=inputFile)
        # A snapshot is needed in order to update the velocities
        snapshot = systemXML.take_snapshot()
        # Assign the required velocities based on the requested temperature
        updatedSnapshot = self.initializeVelocities(snapshot)
        # Finally, restore the snapshot
        systemXML.restore_snapshot(updatedSnapshot)
        return systemXML

    def initializeVelocities(self, snapshot):
        v = np.random.random((len(snapshot.particles.velocity), 3))
        v -= 0.5
        meanv = np.mean(v, 0)
        meanv2 = np.mean(v ** 2, 0)
        fs = np.sqrt(self.temperature / meanv2
        # Shift the velocities such that the average is zero
        v = (v - meanv)
        # Scale the velocities to match the required temperature
        v *= fs
        # Assign the velocities for this MD phase
        snapshot.particles.velocity[:] = v[:]
        return snapshot

    def optimiseStructure(self):
        # Activate the dumping of the trajectory dcd file
        if self.DCDFileWrite is True:
            if self.DCDFileDumpstep is not 0:
                print("Setting DCD dump step to", self.DCDFileDumpstep)
                self.dumpDCD = dump.dcd(filename=self.outputFile.replace('xml', 'dcd'), period=self.DCDFileDumpstep, overwrite=True)
            else:
                self.dumpDCD = dump.dcd(filename=self.outputFile.replace('xml', 'dcd'), period=self.duration / 100.0, overwrite=True)
        else:
            self.dumpDCD = None
        # Set the integrators, groups and timestep
        self.step = integrate.mode_standard(dt=self.timestep)
        self.rigidInt = integrate.nvt_rigid(group=self.rigidGroup, T=self.temperature, tau=self.tau)
        self.nonRigidInt = integrate.nvt(group=self.nonRigidGroup, T=self.temperature, tau=self.tau)
        # Overwrite the log file if this is the first phase, otherwise append to the previous log
        if self.phaseNumber == 0:
            logOverwrite = True
        else:
            logOverwrite = False
        self.energyLog = analyze.log(filename=self.outputLogFileName, quantities=self.logQuantities, period=self.duration / 1000.0, overwrite=logOverwrite)
        callback = None
        # Set up the callback function if the termination condition is not maxt
        if self.terminationCondition == 'KEmin':
            self.loadFromSnapshot = False
            self.lowestKE = 9e999
            self.KEIncreased = 0
            self.firstKEValue = True
            callback = analyze.callback(callback=self.checkKE, period=self.duration / 1000.0)
        print("---=== BEGINNING MOLECULAR DYNAMICS PHASE", self.phaseNumber + 1, "===---")
        # Run the MD simulation
        try:
            run(self.duration)
        except ExitHoomd as exitMessage:
            print(exitMessage)
        # Load the snapshot if required
        if self.terminationCondition == 'KEmin':
            if self.loadFromSnapshot is True:
                print("Loading from snapshot...")
                self.system.restore_snapshot(self.snapshotToLoad)
                del self.snapshotToLoad
        # Create the output XML file
        self.dumpXML = dump.xml(filename=self.outputFile, position=True, image=True, type=True, mass=True, diameter=True, body=True, charge=True, bond=True, angle=True, dihedral=True, improper=True)
        # Clean up all references to this simulation
        self.rigidInt.disable()
        self.nonRigidInt.disable()
        del self.system, self.dumpDCD, self.step, self.rigidInt, self.nonRigidInt, self.rigidGroup, self.nonRigidGroup, callback, self.dumpXML, self.pairClass, self.bondClass, self.angleClass, self.dihedralClass, self.improperClass, self.energyLog
        init.reset()

    def checkKE(self, timestepNumber):
        # Query the current kinetic energy of the system through the energyLog
        currentKE = self.energyLog.query('kinetic_energy')
        # For second and subsequent steps
        if self.firstKEValue is False:
            # Check if the current KE is greater than the minimum so far
            if currentKE >= self.lowestKE:
                if self.KEIncreased == 5:
                    # Found the lowest KE point for at least 5 dumpsteps
                    del self.firstKEValue, self.lowestKE, self.KEIncreased
                    raise ExitHoomd("Lowest Energy Condition Met")
                # Increment a counter that indicates how many times the KE has increased since the minimum
                self.KEIncreased += 1
            else:
                # At at least local KE minimum, so store snapshot
                self.KEIncreased = 0
                self.loadFromSnapshot = True
                self.snapshotToLoad = self.system.take_snapshot(all=True)
                self.lowestKE = currentKE
        else:
            # Skip the first check because the KE fluctuates wildly within the first dump step
            self.firstKEValue = False
        return 0

    def getFFCoeffs(self):
        # First find all of the forcefields specified in the par file
        allFFNames = {}
        for CGSite, directory in self.CGToTemplateDirs.items():
            FFLoc = directory + '/' + self.CGToTemplateForceFields[CGSite]
            if FFLoc not in list(allFFNames.values()):
                allFFNames[CGSite] = FFLoc
        FFList = []
        # Then load in all of the FFs with the appropriate mappings
        for CGSite in list(allFFNames.keys()):
            FFList.append(helperFunctions.loadFFXML(allFFNames[CGSite], mapping = self.newTypeMappings[CGSite]))
        # Combine all of the individual, mapped FFs into one master field
        masterFF = {}
        for FF in FFList:
            for FFType in list(FF.keys()):
                if FFType not in list(masterFF.keys()):
                    masterFF[FFType] = FF[FFType]
                else:
                    masterFF[FFType] += FF[FFType]
        # Finally, assign the expected variables to each value in the masterFF
        self.ljCoeffs = masterFF['lj']
        self.dpdCoeffs = masterFF['dpd']
        self.bondCoeffs = masterFF['bond']
        self.angleCoeffs = masterFF['angle']
        self.dihedralCoeffs = masterFF['dihedral']
        self.improperCoeffs = masterFF['improper']
        # Set Pair Coeffs
        self.pairClass = None
        if self.pairType != 'none':
            # Log the correct pairType energy
            self.logQuantities.append('pair_' + self.pairType + '_energy')
            # HOOMD crashes if you don't specify all pair combinations, so need to make sure we do this.
            atomTypes = sorted(list(set(self.AAMorphologyDict['type'])), key=lambda x: helperFunctions.convertStringToInt(x))
            allPairTypes = []
            # Create a list of all of the pairTypes to ensure that the required coefficients are set
            for atomType1 in atomTypes:
                for atomType2 in atomTypes:
                    pairType = str(atomType1) + '-' + str(atomType2)
                    reversePairType = str(atomType2) + '-' + str(atomType1)
                    if (pairType not in allPairTypes) and (reversePairType not in allPairTypes):
                        allPairTypes.append(pairType)
            # Read in the pairTypes, parameters and coefficients and set them for HOOMD
            if self.pairType == 'dpd':
                self.pairClass = pair.dpd(r_cut=self.pairRCut * self.sScale, T=self.temperature)
                # Use the geometric mixing rule for all possible combinations of the specified forcefield coefficients
                for atomIndex1, atomType1 in enumerate([coeff[0] for coeff in self.dpdCoeffs]):
                    for atomIndex2, atomType2 in enumerate([coeff[0] for coeff in self.dpdCoeffs]):
                        self.pairClass.pair_coeff.set(atomType1, atomType2, A = np.sqrt((self.dpdCoeffs[atomIndex1][1] * self.eScale) * (self.dpdCoeffs[atomIndex2][1] * self.eScale)), r_cut = np.sqrt((self.dpdCoeffs[atomIndex1][2] * self.sScale) * (self.dpdCoeffs[atomIndex2][1] * self.sScale)), gamma = self.pairDPDGammaVal)
                        try:
                            allPairTypes.remove(atomType1 + '-' + atomType2)
                        except:
                            pass
                # Because we've been removing each pair from allPairTypes, all that are left
                # are the pair potentials that are unspecified in the parXX.py (e.g. ghost
                # particle interactions), so set these interactions to zero
                for pairType in allPairTypes:
                    self.pairClass.pair_coeff.set(pairType.split('-')[0], pairType.split('-')[1], A=0.0, r_cut=0.0, gamma=0.0)
            elif self.pairType == 'lj':
                self.pairClass = pair.lj(r_cut=self.pairRCut * self.sScale)
                self.pairClass.set_params(mode='xplor')
                for atomIndex1, atomType1 in enumerate([coeff[0] for coeff in self.ljCoeffs]):
                    for atomIndex2, atomType2 in enumerate([coeff[0] for coeff in self.ljCoeffs]):
                        self.pairClass.pair_coeff.set(atomType1, atomType2, epsilon = np.sqrt((self.ljCoeffs[atomIndex1][1] * self.eScale) * (self.ljCoeffs[atomIndex2][1] * self.eScale)), sigma = np.sqrt((self.ljCoeffs[atomIndex1][2] * self.sScale) * (self.ljCoeffs[atomIndex2][2] * self.sScale)))
                        try:
                            allPairTypes.remove(atomType1 + '-' + atomType2)
                        except:
                            pass
                # Because we've been removing each pair from allPairTypes, all that are left
                # are the pair potentials that are unspecified in the parXX.py (e.g. ghost
                # particle interactions), so set these interactions to zero
                for pairType in allPairTypes:
                    self.pairClass.pair_coeff.set(pairType.split('-')[0], pairType.split('-')[1], epsilon=0.0, sigma=0.0)
            else:
                raise SystemError('Non-DPD/LJ pair potentials not yet hard-coded! Please describe how to interpret them on this line.')
        # Set Bond Coeffs
        # Real bonds
        if self.bondType == 'harmonic':
            self.bondClass = bond.harmonic()
            for bondCoeff in self.bondCoeffs:
                # [k] = kcal mol^{-1} \AA^{-2} * episilon/sigma^{2}, [r0] = \AA * sigma^{2}
                self.bondClass.bond_coeff.set(bondCoeff[0], k=bondCoeff[1] * (self.eScale / (self.sScale**2)), r0=bondCoeff[2] * self.sScale)
        # Ghost bonds
        # If there is no anchoring, rather than change the XML, just set the bond k values to 0.
            if self.groupAnchoring == 'all':
                groupAnchoringTypes = ['X' + CGType for CGType in list(self.CGToTemplateAAIDs.keys())]
            elif self.groupAnchoring == 'none':
                groupAnchoringTypes = []
            else:
                groupAnchoringTypes = ['X' + CGType for CGType in self.groupAnchoring.split(',')]
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
                self.bondClass.bond_coeff.set(bondType, k=1E6, r0=0)
            for bondType in noAnchorBondTypes:
                self.bondClass.bond_coeff.set(bondType, k=0, r0=0)
        else:
            raise SystemError('Non-harmonic bond potentials not yet hard-coded! Please describe how to interpret them on this line.')
        # Set Angle Coeffs
        if self.angleType == 'harmonic':
            self.angleClass = angle.harmonic()
            for angleCoeff in self.angleCoeffs:
                # [k] = kcal mol^{-1} rad^{-2} * epsilon, [t] = rad
                self.angleClass.set_coeff(angleCoeff[0], k=angleCoeff[1] * self.eScale, t0=angleCoeff[2])
        else:
            raise SystemError('Non-harmonic angle potentials not yet hard-coded! Please describe how to interpret them on this line.')
        # Set Dihedral Coeffs
        if self.dihedralType == 'table':
            self.dihedralClass = dihedral.table(width=1000)
            for dihedralCoeff in self.dihedralCoeffs:
                self.dihedralClass.dihedral_coeff.set(dihedralCoeff[0], func=multiHarmonicTorsion, coeff=dict(V0=dihedralCoeff[1] * self.eScale, V1=dihedralCoeff[2] * self.eScale, V2=dihedralCoeff[3] * self.eScale, V3=dihedralCoeff[4] * self.eScale, V4=dihedralCoeff[5] * self.eScale))
        elif self.dihedralType == 'opls':
            self.dihedralClass = dihedral.opls()
            for dihedralCoeff in self.dihedralCoeffs:
                self.dihedralClass.set_coeff(dihedralCoeff[0], k1=dihedralCoeff[1] * self.eScale, k2=dihedralCoeff[2] * self.eScale, k3=dihedralCoeff[3] * self.eScale, k4=dihedralCoeff[4] * self.eScale)
        else:
            raise SystemError('Non-tabulated dihedral potentials not yet hard-coded! Please describe how to interpret them on this line.')
        # Set Improper Coeffs
        self.improperClass = None
        if len(self.improperCoeffs) > 0:
            self.improperClass = improper.harmonic()
            for improperCoeff in self.improperCoeffs:
                self.improperClass.improper_coeff.set(improperCoeff[0], k=improperCoeff[1] * self.eScale, chi=improperCoeff[2])

    def getIntegrationGroups(self):
        # Based on input parameter, return all non-rigid and rigid atoms to be integrated over
        if self.integrationTarget == 'all':
            integrationTypes = list(self.CGToTemplateAAIDs.keys())
        else:
            integrationTypes = self.integrationTarget.split(',')
        # Add in any rigid ghost particles that might need to be integrated too
        ghostIntegrationTypes = ['R' + typeName for typeName in integrationTypes]
        atomIDsToIntegrate = []
        for molecule in self.CGToAAIDMaster:
            for CGSiteID in list(molecule.keys()):
                if molecule[CGSiteID][0] in integrationTypes:
                    atomIDsToIntegrate += molecule[CGSiteID][1]
        for atomID, atomType in enumerate(self.AAMorphologyDict['type']):
            if atomType in ghostIntegrationTypes:
                atomIDsToIntegrate.append(atomID)
        # Create the integrateGroup which contains all of the atoms to be integrated.
        # The rigidGroup constains the intersection of integrateGroup and group.rigid()
        # The nonRigidGroup contains the remainder of atoms in integrateGroup
        integrateGroup = group.tag_list(name="integrateGroup", tags=atomIDsToIntegrate)
        rigidGroup = group.intersection(name="rigidGroup", a=group.rigid(), b=integrateGroup)
        nonRigidGroup = group.difference(name="nonRigidGroup", a=integrateGroup, b=rigidGroup)
        return rigidGroup, nonRigidGroup, integrationTypes


def multiHarmonicTorsion(theta, V0, V1, V2, V3, V4):
    # Definition of multiharmonic dihedral equation based on 5 input parameters to be used by HOOMD
    V = V0 + V1 * np.cos(theta) + V2 * ((np.cos(theta))**2) + V3 * ((np.cos(theta))**3) + V4 * ((np.cos(theta))**4)
    F = V1 * np.sin(theta) + 2 * V2 * np.cos(theta) * np.sin(theta) + 3 * V3 * ((np.cos(theta))**2) * np.sin(theta) + 4 * V4 * ((np.cos(theta))**3) * np.sin(theta)
    return (V, F)


def obtainScaleFactors(parameterDict):
    print("Obtaining correct scaling for epsilon and sigma...")
    # The scaling factors are 1/largestSigma in the LJ coeffs, and 1/largestEpsilon
    LJFFs = []
    for CGSite, directory in parameterDict['CGToTemplateDirs'].items():
        FFLoc = directory + '/' + parameterDict['CGToTemplateForceFields'][CGSite]
        FF = helperFunctions.loadFFXML(FFLoc)
        LJFFs += FF['lj']
    largestSigma = max(list(map(float, np.array(LJFFs)[:, 2])))
    largestEpsilon = max(list(map(float, np.array(LJFFs)[:, 1])))
    return 1 / float(largestSigma), 1 / float(largestEpsilon)


def scaleMorphology(initialMorphology, parameterDict, sScale, eScale):
    # If sScale != 1.0, then scale the morphology and rewrite the phase0 xml
    print("Scaling morphology by sigma =", str(1 / sScale) + "...")
    if sScale != 1.0:
        helperFunctions.scale(initialMorphology, sScale)
    helperFunctions.writeMorphologyXML(initialMorphology, parameterDict['outputMorphDir'] + '/' + parameterDict['morphology'][:-4] + '/morphology/phase0_' + parameterDict['morphology'])


def execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList):
    # Main execution function for runHoomd that performs the required MD phases
    # First, scale the input morphology based on the pair potentials such that
    # the distances and energies are normalised to the strongest pair interaction
    # and the diameter of the largest atom (makes it easier on HOOMDs calculations
    # and ensures that T = 1.0 is an interesting temperature threshold)
    currentFiles = os.listdir(parameterDict['outputMorphDir'] + '/' + parameterDict['morphology'][:-4] + '/morphology')
    # sScale, eScale = obtainScaleFactors(parameterDict)
    print("Under the hood eScaling and sScaling has been disabled.")
    sScale = 1.0
    eScale = 1.0
    # Only scale the morphology if it hasn't been already
    if (parameterDict['overwriteCurrentData'] is False) and ('phase0_' + parameterDict['morphology'] in currentFiles):
        pass
    else:
        scaleMorphology(AAMorphologyDict, parameterDict, sScale, eScale)
    # Reset logfile
    try:
        os.remove(parameterDict['outputMorphDir'] + '/' + parameterDict['morphology'][:-4] + '/morphology/energies_' + parameterDict['morphology'][:-4] + '.log')
    except OSError:
        pass
    # Perform each molecular dynamics phase as specified in the parXX.py
    for phaseNo in range(parameterDict['numberOfPhases']):
        inputFile = 'phase' + str(phaseNo) + '_' + parameterDict['morphology']
        outputFile = 'phase' + str(phaseNo + 1) + '_' + parameterDict['morphology']
        if outputFile in currentFiles:
            if parameterDict['overwriteCurrentData'] is False:
                print(outputFile, "already exists. Skipping...")
                continue
        MDPhase(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, phaseNo, parameterDict['outputMorphDir'] + '/' + parameterDict['morphology'][:-4] + '/morphology/' + inputFile, parameterDict['outputMorphDir'] + '/' + parameterDict['morphology'][:-4] + '/morphology/' + outputFile, sScale, eScale).optimiseStructure()
    finalXMLName = parameterDict['outputMorphDir'] + '/' + parameterDict['morphology'][:-4] + '/morphology/final_' + parameterDict['morphology']
    if 'final_' + parameterDict['morphology'] not in currentFiles:
        # Now all phases are complete, remove the ghost particles from the system
        print("Removing ghost particles to create final output...")
        removeGhostParticles(parameterDict['outputMorphDir'] + '/' + parameterDict['morphology'][:-4] + '/morphology/' + outputFile, finalXMLName, sigma = sScale)
    # Finally, update the pickle file with the most recent and realistic
    # AAMorphologyDict so that we can load it again further along the pipeline
    AAMorphologyDict = helperFunctions.loadMorphologyXML(finalXMLName)
    # Now that we've obtained the final fine-grained morphology, we need to fix the images to prevent
    # issues with obtaining the chromophores and running them through the ZINDO/S calculations later...
    AAMorphologyDict = helperFunctions.fixImages(AAMorphologyDict)
    # ...add in the unwrapped positions...
    AAMorphologyDict = helperFunctions.addUnwrappedPositions(AAMorphologyDict)
    # ...and write the pickle file.
    helperFunctions.writePickle((AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList), parameterDict['outputMorphDir'] + '/' + parameterDict['morphology'][:-4] + '/code/' + parameterDict['morphology'][:-4] + '.pickle')
    return AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList


def removeGhostParticles(lastPhaseXML, outputFileName, sigma = 1.0):
    # Remove all the ghost particles from the morphology for the final output
    finalMorphology = helperFunctions.loadMorphologyXML(lastPhaseXML)
    # Determine the atomIDs for each particle beginning with the letters 'X'
    # or 'R' - these are the ghost particles
    atomIDsToRemove = []
    for atomID, atomType in enumerate(finalMorphology['type']):
        if (atomType[0] == 'X') or (atomType[0] == 'R'):
            # This is a ghost particle
            atomIDsToRemove.append(atomID)
    # Reverse sorting trick so that the location indices don't change as we
    # delete particles from the system
    atomIDsToRemove.sort(reverse=True)
    # Now delete the atoms from the morphology
    atomAttribs = ['position', 'image', 'type', 'mass', 'diameter', 'body', 'charge']
    for atomID in atomIDsToRemove:
        for key in atomAttribs:
            finalMorphology[key].pop(atomID)
    finalMorphology['natoms'] -= len(atomIDsToRemove)
    # Delete any constraints associated with those atoms that have been removed
    atomConstraints = ['bond', 'angle', 'dihedral', 'improper']
    for key in atomConstraints:
        constraintsToRemove = []
        for constraintNo, constraint in enumerate(finalMorphology[key]):
            for atomID in constraint[1:]:
                if (atomID in atomIDsToRemove) and (constraintNo not in constraintsToRemove):
                    constraintsToRemove.append(constraintNo)
        constraintsToRemove.sort(reverse=True)
        for constraintNo in constraintsToRemove:
            finalMorphology[key].pop(constraintNo)
    # Output the final morphology
    helperFunctions.writeMorphologyXML(finalMorphology, outputFileName, sigma)


if __name__ == "__main__":
    try:
        pickleFile = sys.argv[1]
    except:
        print("Please specify the pickle file to load to continue the pipeline from this point.")
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = helperFunctions.loadPickle(pickleFile)
    execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList)
