from hoomd_script import *
import numpy as np
import copy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cPickle as pickle
import helperFunctions
import sys

def multiHarmonicTorsion(theta, V0, V1, V2, V3, V4):
    V = V0 + V1*np.cos(theta) + V2*((np.cos(theta))**2) + V3*((np.cos(theta))**3) + V4*((np.cos(theta))**4)
    F = V1*np.sin(theta) + 2*V2*np.cos(theta)*np.sin(theta) + 3*V3*((np.cos(theta))**2)*np.sin(theta) + 4*V4*((np.cos(theta))**3)*np.sin(theta)
    return (V, F)


class ExitHoomd(Exception):
    def __init__(self, string, moleculeName):
        self.string = string+" At Timestep = "+str(get_step())+" For Molecule = "+moleculeName
    def __str__(self):
        return self.string


class hoomdRun:
    def __init__(self, fileName, CGMoleculeDict, CG2AAIDs, eScale, sScale, continueData, AAMorphologyDict):
        slashLocs = helperFunctions.findIndex(fileName, '/')
        underscoreLocs = helperFunctions.findIndex(fileName, '_')
        self.morphologyName = fileName[underscoreLocs[-1]+1:-4]
        self.saveDirectory = fileName[:slashLocs[-1]]+'/'
        self.fileName = fileName
        self.CGMoleculeDict = CGMoleculeDict
        self.AAMorphologyDict = AAMorphologyDict
        self.CG2AAIDs = CG2AAIDs
        self.eScale = eScale
        self.sScale = sScale
        self.dumpPeriod = 1e2
        self.mainTrajDumpPeriod = 1e2
        self.T = 1.0
        self.tau = 1.0
        #self.dtPhase1 = 1e-5

        # Phase 1 = KE-truncated relaxation with pair potentials off
        self.dtPhase1 = 1e-3 # No Pairs
        self.dtPhase2 = 1e-3 # DPD Hydrogens and Sidechains, gradientRamp = 100, A = 1*eScale for all C*-C* and H1-*, 0 for else. r_cut = 2.96*sScale, gamma = 0
        self.dtPhase3 = 1e-9 # LJ
        self.dtPhase4 = 1e-7 # LJ
        self.dtPhase5 = 1e-6 # LJ
        self.dtPhase6 = 1e-5 # LJ
        self.phase1RunLength = 1e5 # Maximum, this run is KE truncated
        self.phase2RunLength = 1e4
        self.phase3RunLength = 1e2 
        self.phase4RunLength = 1e2
        self.phase5RunLength = 1e5
        self.phase6RunLength = 1e5
        self.outputXML = self.saveDirectory+'relaxed_'+self.morphologyName+'.xml'
        self.outputDCD = self.saveDirectory+'relaxed_'+self.morphologyName+'.dcd'
        self.outputLOG = self.saveDirectory+'energies_'+self.morphologyName+'.log'
        self.eScale = eScale
        self.sScale = sScale
        self.overwriteEnergies = True
        # self.previousTempXML = False
        # Check the save directory for previous runs so that we can continue where we left off
        self.runPhase1 = continueData[0]
        self.runPhase2 = continueData[1]
        self.runPhase3 = continueData[2]
        self.runPhase4 = continueData[3]
        self.runPhase5 = continueData[4]
        self.runPhase6 = continueData[5]
        self.continuePhase6 = continueData[6]
        self.continueFile = continueData[7]
        

    def initialiseRun(self, inputFilename, pairType='lj', gradientRamp=1.0, rigidBodies=True):
        if rigidBodies == True:
            # Sort out the groups. Currently self.thioGroupIDs contains all of the
            # Hydrogens in the thiophene ring and at the end of the molecules too.
            # We don't want these to be rigid, so want them to be added to one of
            # the sidechain groups self.alk1Group or self.alk2Group instead.
            self.system = init.read_xml(filename=inputFilename)
            hydrogens = group.type(name="hydrogens", type='H1')
            fullThiophenes = group.tag_list(name="thioIncH", tags=self.thioGroupIDs)
            self.thioGroup = group.difference(name="thio", a = fullThiophenes, b = hydrogens)
            alk1GroupWithoutH = group.tag_list(name="alk1WOH", tags=self.alk1GroupIDs)
            alk1Group = group.union(name="alk1IncExtraH", a = alk1GroupWithoutH, b = hydrogens)
            alk2Group = group.tag_list(name="alk2", tags=self.alk2GroupIDs)
            self.sideChainsGroup = group.union(name="sideChains", a = alk1Group, b = alk2Group)
        else:
            # Current XML has all the rigid body information in it. Remove this data.
            morphologyWithRigidBodies = helperFunctions.loadMorphologyXML(inputFilename)
            morphologyWithoutRigidBodies = helperFunctions.removeRigidBodies(morphologyWithRigidBodies)
            helperFunctions.writeMorphologyXML(morphologyWithoutRigidBodies, inputFilename.replace('phase1', 'phase1Flex'))
            self.system = init.read_xml(filename=inputFilename.replace('phase1', 'phase1Flex'))
            self.thioGroup = group.tag_list(name="thio", tags=self.thioGroupIDs)
            self.alk1Group = group.tag_list(name="alk1", tags=self.alk1GroupIDs)
            self.alk2Group = group.tag_list(name="alk2", tags=self.alk2GroupIDs)
        self.setForcefieldCoeffs(self.eScale, self.sScale, pairType, gradientRamp)
        if pairType == 'none':
            self.energyLog = analyze.log(filename = self.outputLOG, quantities=['potential_energy', 'kinetic_energy', 'bond_harmonic_energy', 'angle_harmonic_energy', 'dihedral_table_energy', 'temperature', 'pressure', 'volume'], period=self.dumpPeriod, overwrite = self.overwriteEnergies)
        elif pairType == 'dpd':
            self.energyLog = analyze.log(filename = self.outputLOG, quantities=['potential_energy', 'kinetic_energy', 'pair_dpd_energy', 'bond_harmonic_energy', 'angle_harmonic_energy', 'dihedral_table_energy', 'temperature', 'pressure', 'volume'], period=self.dumpPeriod, overwrite = self.overwriteEnergies)
        elif pairType == 'lj':
            self.energyLog = analyze.log(filename = self.outputLOG, quantities=['potential_energy', 'kinetic_energy', 'pair_lj_energy', 'bond_harmonic_energy', 'angle_harmonic_energy', 'dihedral_table_energy', 'temperature', 'pressure', 'volume'], period=self.dumpPeriod, overwrite = self.overwriteEnergies)
        self.overwriteEnergies = False


    def setForcefieldCoeffs(self, eScale, sScale, pairType, gradientRamp):
        # Forcefield parameters obtained from:
        # Bhatta, R. S., Yimer, Y. Y, Perry, D. S., Tsige, M., "Improved Force Field for Molecular Modeling of Poly(3-Hexylthiophene)", 2013, J. Phys. Chem. B, DOI: 10.1021/jp404629a
        # Marcon, V., Raos, G., "Free Energies of Molecular Crystal Surfaces by Computer Simlations: Application to Tetrathiophene", 2006, J. Am. Chem. Soc., DOI: 10.1021/ja056548t
        # Marcon, V., Raos, G., "Molecular Modeling of Crystalline Oligothiophenes: Testing and Development of Improved Force Fields", 2004, J. Phys. Chem. B, DOI: 10.1021/jp047128d
        # Huang, D. M., Faller, R., Do, K., Moule, A. J., "Coarse-Grained Computer Simulations of Polymer/Fullerene Bulk Heterojunctions for Organic Photovoltaic Applications", 2010, J. Chem. Theory Comput., DOI: 10.1021/ct900496t
        # Jorgensen, W. L., Maxwell, D. S., Tirdao-Rives, J., "Development and Testing of the OPLS All-Atom Forcefield on Conformational Energetics and Properties of Organic Liquids", 1996, J. Am. Chem. Soc. DOI: 10.1021/ja9621760
        # Bhatta, R. S., Yimer, Y. Y., Tsige, M., Perry, D. S., "Conformations and Torsional Potentials of Poly(3-Hexylthiophene) Oligomers: Density Functional Calculations Up to the Dodecamer", 2012, Comput. & Theor. Chem., DOI: 10.1016/j.comptc.2012.06.026
        if pairType == 'none':
            self.setLJParameters(eScale, sScale, 0)
        elif pairType == 'dpd':
            self.setDPDParameters(eScale, sScale, gradientRamp)
        elif pairType == 'lj':
            self.setLJParameters(eScale, sScale, gradientRamp)
        self.setHarmonicBondParameters(eScale, sScale)
        self.setHarmonicAngleParameters(eScale, sScale)
        self.setTableDihedralParameters(eScale, sScale)
        self.setImproperParameters(eScale, sScale)

        
    def setDPDParameters(self, eScale, sScale, gradientRamp):
        gammaVal = 0.0
        eScale *= gradientRamp
        self.pair = pair.dpd(r_cut=10*sScale, T=1.0)
        sScale *= 2**(1./6.)
        for type1 in ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'S1', 'H1', 'T', 'X1', 'X2', 'X3']:
            for type2 in ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'S1', 'H1', 'T', 'X1', 'X2', 'X3']:
                if ('C' in type1) and ('C' in type2):
                    self.pair.pair_coeff.set(type1, type2, A = 1*eScale, r_cut = 2.96*sScale, gamma=gammaVal)
                else:
                    self.pair.pair_coeff.set(type1, type2, A = 0, r_cut = 0, gamma=gammaVal)
        for type1 in ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'S1', 'H1', 'T', 'X1', 'X2', 'X3']:
            self.pair.pair_coeff.set(type1, 'H1', A = 1*eScale, r_cut = 2.96*sScale, gamma=gammaVal)
        self.pair.pair_coeff.set('H1', 'S1', A = 1*eScale, r_cut = 2.96*sScale, gamma=gammaVal)
        # self.pair.pair_coeff.set('C1','H1',A=0.046*eScale,r_cut=2.979*sScale,gamma=gammaVal)
        # self.pair.pair_coeff.set('C2','H1',A=0.046*eScale,r_cut=2.979*sScale,gamma=gammaVal)
        # self.pair.pair_coeff.set('C3','H1',A=0.044*eScale,r_cut=2.958*sScale,gamma=gammaVal)
        # self.pair.pair_coeff.set('C4','H1',A=0.044*eScale,r_cut=2.958*sScale,gamma=gammaVal)
        # self.pair.pair_coeff.set('C5','H1',A=0.044*eScale,r_cut=2.958*sScale,gamma=gammaVal)
        # self.pair.pair_coeff.set('C6','H1',A=0.044*eScale,r_cut=2.958*sScale,gamma=gammaVal)
        # self.pair.pair_coeff.set('C7','H1',A=0.044*eScale,r_cut=2.958*sScale,gamma=gammaVal)
        # self.pair.pair_coeff.set('C8','H1',A=0.044*eScale,r_cut=2.958*sScale,gamma=gammaVal)
        # self.pair.pair_coeff.set('C9','H1',A=0.046*eScale,r_cut=2.931*sScale,gamma=gammaVal)
        # self.pair.pair_coeff.set('C10','H1',A=0.046*eScale,r_cut=2.979*sScale,gamma=gammaVal)
        # self.pair.pair_coeff.set('H1','H1',A=0.030*eScale,r_cut=2.500*sScale,gamma=gammaVal)
        # self.pair.pair_coeff.set('H1','S1',A=0.087*eScale,r_cut=2.979*sScale,gamma=gammaVal)
        
        
    def setLJParameters(self, eScale, sScale, gradientRamp):
        eScale *= gradientRamp
        self.pair = pair.lj(r_cut=10*sScale)
        self.pair.set_params(mode='xplor')
        self.pair.pair_coeff.set('C1','C1',epsilon=0.070*eScale,sigma=3.550*sScale)
        self.pair.pair_coeff.set('C1','C2',epsilon=0.070*eScale,sigma=3.550*sScale)
        self.pair.pair_coeff.set('C1','C3',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C1','C4',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C1','C5',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C1','C6',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C1','C7',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C1','C8',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C1','C9',epsilon=0.070*eScale,sigma=3.550*sScale)
        self.pair.pair_coeff.set('C1','C10',epsilon=0.070*eScale,sigma=3.550*sScale)
        self.pair.pair_coeff.set('C1','H1',epsilon=0.046*eScale,sigma=2.979*sScale)
        self.pair.pair_coeff.set('C1','S1',epsilon=0.132*eScale,sigma=3.550*sScale)
        self.pair.pair_coeff.set('C2','C2',epsilon=0.070*eScale,sigma=3.550*sScale)
        self.pair.pair_coeff.set('C2','C3',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C2','C4',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C2','C5',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C2','C6',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C2','C7',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C2','C8',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C2','C9',epsilon=0.070*eScale,sigma=3.550*sScale)
        self.pair.pair_coeff.set('C2','C10',epsilon=0.070*eScale,sigma=3.550*sScale)
        self.pair.pair_coeff.set('C2','H1',epsilon=0.046*eScale,sigma=2.979*sScale)
        self.pair.pair_coeff.set('C2','S1',epsilon=0.132*eScale,sigma=3.550*sScale)
        self.pair.pair_coeff.set('C3','C3',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.pair.pair_coeff.set('C3','C4',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.pair.pair_coeff.set('C3','C5',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.pair.pair_coeff.set('C3','C6',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.pair.pair_coeff.set('C3','C7',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.pair.pair_coeff.set('C3','C8',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.pair.pair_coeff.set('C3','C9',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C3','C10',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C3','H1',epsilon=0.044*eScale,sigma=2.958*sScale)
        self.pair.pair_coeff.set('C3','S1',epsilon=0.128*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C4','C4',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.pair.pair_coeff.set('C4','C5',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.pair.pair_coeff.set('C4','C6',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.pair.pair_coeff.set('C4','C7',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.pair.pair_coeff.set('C4','C8',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.pair.pair_coeff.set('C4','C9',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C4','C10',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C4','H1',epsilon=0.044*eScale,sigma=2.958*sScale)
        self.pair.pair_coeff.set('C4','S1',epsilon=0.128*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C5','C5',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.pair.pair_coeff.set('C5','C6',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.pair.pair_coeff.set('C5','C7',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.pair.pair_coeff.set('C5','C8',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.pair.pair_coeff.set('C5','C9',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C5','C10',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C5','H1',epsilon=0.044*eScale,sigma=2.958*sScale)
        self.pair.pair_coeff.set('C5','S1',epsilon=0.128*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C6','C6',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.pair.pair_coeff.set('C6','C7',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.pair.pair_coeff.set('C6','C8',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.pair.pair_coeff.set('C6','C9',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C6','C10',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C6','H1',epsilon=0.044*eScale,sigma=2.958*sScale)
        self.pair.pair_coeff.set('C6','S1',epsilon=0.128*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C7','C7',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.pair.pair_coeff.set('C7','C8',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.pair.pair_coeff.set('C7','C9',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C7','C10',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C7','H1',epsilon=0.044*eScale,sigma=2.958*sScale)
        self.pair.pair_coeff.set('C7','S1',epsilon=0.128*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C8','C8',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.pair.pair_coeff.set('C8','C9',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C8','C10',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C8','H1',epsilon=0.044*eScale,sigma=2.958*sScale)
        self.pair.pair_coeff.set('C8','S1',epsilon=0.128*eScale,sigma=3.525*sScale)
        self.pair.pair_coeff.set('C9','C9',epsilon=0.070*eScale,sigma=3.550*sScale)
        self.pair.pair_coeff.set('C9','C10',epsilon=0.070*eScale,sigma=3.550*sScale)
        self.pair.pair_coeff.set('C9','H1',epsilon=0.046*eScale,sigma=2.931*sScale)
        self.pair.pair_coeff.set('C9','S1',epsilon=0.132*eScale,sigma=3.550*sScale)
        self.pair.pair_coeff.set('C10','C10',epsilon=0.070*eScale,sigma=3.550*sScale)
        self.pair.pair_coeff.set('C10','H1',epsilon=0.046*eScale,sigma=2.979*sScale)
        self.pair.pair_coeff.set('C10','S1',epsilon=0.132*eScale,sigma=3.550*sScale)
        self.pair.pair_coeff.set('H1','H1',epsilon=0.030*eScale,sigma=2.500*sScale)
        self.pair.pair_coeff.set('H1','S1',epsilon=0.087*eScale,sigma=2.979*sScale)
        self.pair.pair_coeff.set('S1','S1',epsilon=0.250*eScale,sigma=3.550*sScale)
        # Set Ghost Particle Interactions
        for ghostAtomType in ['T', 'X1', 'X2', 'X3']:
            for realAtomType in ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'S1', 'H1', 'T', 'X1', 'X2', 'X3']:
                self.pair.pair_coeff.set(ghostAtomType, realAtomType, epsilon=0, sigma = 0)
        print self.eScale, eScale

        
    def setHarmonicBondParameters(self, eScale, sScale):
        self.b = bond.harmonic()
        # Ring Bonds [k] = kcal mol^{-1} \AA^{-2} * episilon/sigma^{2}, [r] = \AA * sigma^{2}
        self.b.bond_coeff.set('C1-S1',k=582.50*(eScale/(sScale**2)),r0=1.73373*sScale)
        self.b.bond_coeff.set('C10-S1',k=582.50*(eScale/(sScale**2)),r0=1.73373*sScale)
        self.b.bond_coeff.set('C1-C2',k=1028.54*(eScale/(sScale**2)),r0=1.37368*sScale) # THIS BOND IS A DOUBLE
        self.b.bond_coeff.set('C2-C9',k=906.2*(eScale/(sScale**2)),r0=1.43277*sScale)
        self.b.bond_coeff.set('C9-H1',k=741.26*(eScale/(sScale**2)),r0=1.0822*sScale)
        self.b.bond_coeff.set('C10-C9',k=1028.54*(eScale/(sScale**2)),r0=1.37368*sScale) # THIS BOND IS A DOUBLE
        # Alkyl Bonds
        self.b.bond_coeff.set('C2-C3',k=599.64*(eScale/(sScale**2)),r0=1.50884*sScale)
        self.b.bond_coeff.set('C3-H1',k=655.09*(eScale/(sScale**2)),r0=1.09827*sScale)
        self.b.bond_coeff.set('C3-C4',k=536.00*(eScale/(sScale**2)),r0=1.54158*sScale)
        self.b.bond_coeff.set('C4-H1',k=680.00*(eScale/(sScale**2)),r0=1.09527*sScale)
        self.b.bond_coeff.set('C4-C5',k=536.00*(eScale/(sScale**2)),r0=1.54158*sScale)
        self.b.bond_coeff.set('C5-H1',k=680.00*(eScale/(sScale**2)),r0=1.09527*sScale)
        self.b.bond_coeff.set('C5-C6',k=536.00*(eScale/(sScale**2)),r0=1.54158*sScale)
        self.b.bond_coeff.set('C6-H1',k=680.00*(eScale/(sScale**2)),r0=1.09527*sScale)
        self.b.bond_coeff.set('C6-C7',k=536.00*(eScale/(sScale**2)),r0=1.54158*sScale)
        self.b.bond_coeff.set('C7-H1',k=680.00*(eScale/(sScale**2)),r0=1.09527*sScale)
        self.b.bond_coeff.set('C7-C8',k=536.00*(eScale/(sScale**2)),r0=1.54158*sScale)
        self.b.bond_coeff.set('C8-H1',k=680.00*(eScale/(sScale**2)),r0=1.09527*sScale)
        # Inter-monomer Bonds
        self.b.bond_coeff.set('C1-C10',k=784.58*(eScale/(sScale**2)),r0=1.45*sScale)
        # Terminating Bonds
        self.b.bond_coeff.set('C1-H1',k=741.26*(eScale/(sScale**2)),r0=1.0822*sScale)
        self.b.bond_coeff.set('C10-H1',k=741.26*(eScale/(sScale**2)),r0=1.0822*sScale)
        # Ghost-Anchor Bonds
        self.b.bond_coeff.set('T-X1',k=1E6,r0=0)
        self.b.bond_coeff.set('X2-C4',k=1E6,r0=0)
        self.b.bond_coeff.set('X3-C7',k=1E6,r0=0)

        
    def setHarmonicAngleParameters(self, eScale, sScale):
        self.a = angle.harmonic()
        # Ring Bond Angles [k] = kcal mol^{-1} rad^{-2} * epsilon, [t] = rad
        self.a.set_coeff('C10-C9-C2',k=39.582*eScale,t0=1.97784)
        self.a.set_coeff('C10-C9-H1',k=35.263*eScale,t0=2.14639)
        self.a.set_coeff('C10-S1-C1',k=86.36*eScale,t0=1.61921)
        self.a.set_coeff('C9-C10-S1',k=86.36*eScale,t0=1.92496)
        self.a.set_coeff('C1-C2-C9',k=39.582*eScale,t0=1.97784)
        self.a.set_coeff('C2-C9-H1',k=35.263*eScale,t0=2.15897)
        self.a.set_coeff('C2-C1-S1',k=86.36*eScale,t0=1.92496)
        # Alkyl Bond Angles
        self.a.set_coeff('C3-C2-C9',k=166.545*eScale,t0=2.15335)
        self.a.set_coeff('C2-C3-C4',k=120.14*eScale,t0=2.01481)

        self.a.set_coeff('C1-C2-C3',k=166.32*eScale,t0=2.17388)
        self.a.set_coeff('C3-C4-C5',k=58.35*eScale,t0=1.96699)
        self.a.set_coeff('C4-C5-C6',k=58.35*eScale,t0=1.96699)
        self.a.set_coeff('C5-C6-C7',k=58.35*eScale,t0=1.96699)
        self.a.set_coeff('C6-C7-C8',k=58.35*eScale,t0=1.96699)
        # Alkyl Bond Angles, with Hydrogens Only
        self.a.set_coeff('C2-C3-H1',k=74.06*eScale,t0=1.90571)
        self.a.set_coeff('C3-C4-H1',k=37.5*eScale,t0=1.93208)
        self.a.set_coeff('C4-C3-H1',k=37.5*eScale,t0=1.93208)
        self.a.set_coeff('C4-C5-H1',k=37.5*eScale,t0=1.93208)
        self.a.set_coeff('C5-C4-H1',k=37.5*eScale,t0=1.93208)
        self.a.set_coeff('C5-C6-H1',k=37.5*eScale,t0=1.93208)
        self.a.set_coeff('C6-C5-H1',k=37.5*eScale,t0=1.93208)
        self.a.set_coeff('C6-C7-H1',k=37.5*eScale,t0=1.93208)
        self.a.set_coeff('C7-C6-H1',k=37.5*eScale,t0=1.93208)
        self.a.set_coeff('C7-C8-H1',k=37.5*eScale,t0=1.93208)
        self.a.set_coeff('C8-C7-H1',k=37.5*eScale,t0=1.93208)
        self.a.set_coeff('H1-C8-H1',k=33.0*eScale,t0=1.88146)
        self.a.set_coeff('H1-C7-H1',k=33.0*eScale,t0=1.88146)
        self.a.set_coeff('H1-C6-H1',k=33.0*eScale,t0=1.88146)
        self.a.set_coeff('H1-C5-H1',k=33.0*eScale,t0=1.88146)
        self.a.set_coeff('H1-C4-H1',k=33.0*eScale,t0=1.88146)
        self.a.set_coeff('H1-C3-H1',k=33.0*eScale,t0=1.88146)   
        # Inter-monomer Bond Angles
        self.a.set_coeff('C1-C10-C9',k=54.694*eScale,t0=2.27137)
        self.a.set_coeff('C1-C10-S1',k=41.74*eScale,t0=2.08687)
        self.a.set_coeff('S1-C1-C10',k=41.74*eScale,t0=2.08687)
        self.a.set_coeff('C2-C1-C10',k=54.694*eScale,t0=2.27137)

        
    def setTableDihedralParameters(self, eScale, sScale):
        self.d = dihedral.table(width=1000)
        # Ring Dihedrals
        self.d.dihedral_coeff.set('C10-C2-C9-C1',func=multiHarmonicTorsion, coeff=dict(V0=126.32*eScale,V1=-109.81*eScale,V2=-19.738*eScale,V3=-25.303*eScale,V4=28.53*eScale))
        self.d.dihedral_coeff.set('C10-S1-C1-C2',func=multiHarmonicTorsion, coeff=dict(V0=126.32*eScale,V1=-109.81*eScale,V2=-19.738*eScale,V3=-25.303*eScale,V4=28.53*eScale))
        self.d.dihedral_coeff.set('C1-S1-C10-C9',func=multiHarmonicTorsion, coeff=dict(V0=126.32*eScale,V1=-109.81*eScale,V2=-19.738*eScale,V3=-25.303*eScale,V4=28.53*eScale))
        #self.d.dihedral_coeff.set('C9-C2-C1-S1',func=multiHarmonicTorsion, coeff=dict(V0=126.32,V1=-109.81,V2=-19.738,V3=-25.303,V4=28.53)) # Guessed this one
        self.d.dihedral_coeff.set('C2-C9-C10-S1',func=multiHarmonicTorsion, coeff=dict(V0=126.32*eScale,V1=-109.81*eScale,V2=-19.738*eScale,V3=-25.303*eScale,V4=28.53*eScale))
        # Alkyl Dihedrals
        self.d.dihedral_coeff.set('C10-C9-C2-C3',func=multiHarmonicTorsion, coeff=dict(V0=117.65*eScale,V1=238.26*eScale,V2=205.96*eScale,V3=112.81*eScale,V4=27.467*eScale))
        self.d.dihedral_coeff.set('C4-C3-C2-C9',func=multiHarmonicTorsion, coeff=dict(V0=0.3175*eScale,V1=1.127*eScale,V2=14.143*eScale,V3=-22.297*eScale,V4=6.7188*eScale))
        #self.d.dihedral_coeff.set('C9-C2-C3-H1',func=multiHarmonicTorsion, coeff=dict(V0=117.65,V1=238.26,V2=205.96,V3=112.81,V4=27.467)) # Guessed this one
        self.d.dihedral_coeff.set('C2-C3-C4-C5',func=multiHarmonicTorsion, coeff=dict(V0=2.4469*eScale,V1=-6.3946*eScale,V2=10.747*eScale,V3=30.695*eScale,V4=11.139*eScale))
        self.d.dihedral_coeff.set('C3-C4-C5-C6',func=multiHarmonicTorsion, coeff=dict(V0=1.9475*eScale,V1=-3.7121*eScale,V2=1.388*eScale,V3=8.6305*eScale,V4=1.6008*eScale))
        self.d.dihedral_coeff.set('C4-C5-C6-C7',func=multiHarmonicTorsion, coeff=dict(V0=1.8922*eScale,V1=-3.4904*eScale,V2=1.4665*eScale,V3=7.1418*eScale,V4=0.2859*eScale))
        self.d.dihedral_coeff.set('C5-C6-C7-C8',func=multiHarmonicTorsion, coeff=dict(V0=1.9788*eScale,V1=-3.8476*eScale,V2=1.1614*eScale,V3=7.419*eScale,V4=0.4146*eScale))
        # Inter-monomer Dihedrals
        self.d.dihedral_coeff.set('C1-C10-C9-C2',func=multiHarmonicTorsion, coeff=dict(V0=75.595*eScale,V1=116.*eScale,V2=42.679*eScale,V3=-1.528*eScale,V4=-3.8137*eScale))
        self.d.dihedral_coeff.set('C1-C10-S1-C1',func=multiHarmonicTorsion, coeff=dict(V0=158.7*eScale,V1=418.34*eScale,V2=521.33*eScale,V3=376.73*eScale,V4=115.12*eScale)) # 3, 25, 29, 28
        self.d.dihedral_coeff.set('S1-C1-C10-S1',func=multiHarmonicTorsion, coeff=dict(V0=2.9533*eScale,V1=0.1571*eScale,V2=-4.2326*eScale,V3=0.3979*eScale,V4=1.8855*eScale))
        self.d.dihedral_coeff.set('C2-C1-C10-S1',func=multiHarmonicTorsion, coeff=dict(V0=2.9533*eScale,V1=-0.1571*eScale,V2=-4.2326*eScale,V3=-0.3979*eScale,V4=1.8855*eScale))

        
    def setImproperParameters(self, eScale, sScale):
        # Placeholder for impropers if we use them later
        self.i = None

        
    def optimiseStructure(self):
        thioAtoms = []
        alk1Atoms = []
        alk2Atoms = []

        self.thioGroupIDs = []
        self.alk1GroupIDs = []
        self.alk2GroupIDs = []
        for moleculeNo in range(len(self.CG2AAIDs)):
            for CGAtomID in self.CG2AAIDs[moleculeNo].keys():
                COMPosn = self.CGMoleculeDict['position'][CGAtomID]
                self.CG2AAIDs[moleculeNo][CGAtomID].append(COMPosn)
                # CG2AAIDs is now sorted by the CG Atom ID, and is characterised by:
                # CG2AAIDs[moleculeNumber][CGAtomID] = [CGType, [atom1InGroup, atom2InGroup...], targetCOMPosn]
                # Build the groups:
                if self.CG2AAIDs[moleculeNo][CGAtomID][0] == 'thio':
                    thioAtoms.append(self.CG2AAIDs[moleculeNo][CGAtomID][1:])
                    self.thioGroupIDs += self.CG2AAIDs[moleculeNo][CGAtomID][1]
                elif self.CG2AAIDs[moleculeNo][CGAtomID][0] == 'alk1':
                    alk1Atoms.append(self.CG2AAIDs[moleculeNo][CGAtomID][1:])
                    self.alk1GroupIDs += self.CG2AAIDs[moleculeNo][CGAtomID][1]
                elif self.CG2AAIDs[moleculeNo][CGAtomID][0] == 'alk2':
                    alk2Atoms.append(self.CG2AAIDs[moleculeNo][CGAtomID][1:])
                    self.alk2GroupIDs += self.CG2AAIDs[moleculeNo][CGAtomID][1]
        # Add the 'T' type ghost particle into self.thioGroupIDs
        for atomID, atomType in enumerate(self.AAMorphologyDict['type']):
            if atomType == 'T':
                self.thioGroupIDs.append(atomID)
                    
        if self.runPhase1 == True:
            self.initialiseRun(self.fileName, pairType='none')
            phase1DumpDCD = dump.dcd(filename=self.outputDCD.replace('relaxed', 'phase1'), period=10, overwrite=True)
            phase1Step = integrate.mode_standard(dt = self.dtPhase1)
            # phase1 = integrate.brownian(group=group.all(), seed=3, dscale=1e11, T=self.T)
            #phase1Flex = integrate.nve(group=self.sideChainsGroup, limit=0.001)
            phase1Flex = integrate.nvt(group=self.sideChainsGroup, T=self.T, tau=self.tau)
            phase1Rig = integrate.nvt_rigid(group=self.thioGroup, T=self.T, tau=self.tau)
            self.initialPotentialEnergies = []
            self.initialKineticEnergies = []
            self.loadFromSnapshot = False
            self.lowestKE = 9e999
            self.KEIncreased = 0
            self.firstKEValue = True
            # Modeler_hoomd wipes out the velocity data, and so we start with 0 kinetic energy. Skip this step (firstKEValue == True). Then, take the snapshot with the lowest KE.
            checkKEs = analyze.callback(callback = self.checkKE, period=self.dumpPeriod)
            try:
                run_upto(self.phase1RunLength)
                #run(self.maximumInitialRunLength)
            except ExitHoomd as exitMessage:
                print exitMessage
            if self.loadFromSnapshot == True:
                print "Loading from snapshot..."
                self.system.restore_snapshot(self.snapshotToLoad)
            phase1DumpXML = dump.xml(filename=self.outputXML.replace('relaxed', 'phase1'), position = True, image = True, type = True, mass = True, diameter = True, body = True, charge = True, bond = True, angle = True, dihedral = True, improper = True)
            phase1Flex.disable()
            phase1Rig.disable()
            phase1DumpDCD.disable()
            del self.system, self.thioGroup, self.sideChainsGroup, self.energyLog, self.pair, self.b, self.a, self.d, self.i, phase1DumpDCD, phase1Step, phase1DumpXML, phase1Rig, phase1Flex
            init.reset()
        else:
            print "Phase 1 already completed for this morphology...skipping"

        if self.runPhase2 == True:
            self.initialiseRun(self.outputXML.replace('relaxed', 'phase1'), pairType='dpd', rigidBodies=True, gradientRamp = 1e2)
            phase2DumpDCD = dump.dcd(filename=self.outputDCD.replace('relaxed', 'phase2'), period=100, overwrite=True)
            phase2Step = integrate.mode_standard(dt = self.dtPhase2)
            phase2Hydrogens = integrate.nvt(group=group.type(name="hydrogens", type='H1'), T=self.T, tau=self.tau)
            phase2Sidechains = integrate.nvt(group=group.difference(name="SidechainsWOH", a = self.sideChainsGroup, b = group.type(name="hydrogens", type='H1')), T=self.T, tau=self.tau)
            #phase2Rig = integrate.nvt_rigid(group=self.thioGroup, T=self.T, tau=self.tau)
            run(self.phase2RunLength)
            phase2DumpXML = dump.xml(filename=self.outputXML.replace('relaxed', 'phase2'), position = True, image = True, type = True, mass = True, diameter = True, body = True, charge = True, bond = True, angle = True, dihedral = True, improper = True)
            phase2Hydrogens.disable()
            phase2DumpDCD.disable()
            del self.system, self.thioGroup, self.sideChainsGroup, self.energyLog, self.pair, self.b, self.a, self.d, self.i, phase2DumpDCD, phase2Step, phase2DumpXML, phase2Hydrogens, phase2Sidechains
            init.reset()
        else:
            print "Phase 2 already completed for this morphology...skipping"
        # initDump.disable()
        # debugDump = dump.dcd(filename=self.outputDCD, period=1, overwrite=False)


        if self.runPhase3 == True:
            self.initialiseRun(self.outputXML.replace('relaxed', 'phase2'), pairType='lj', rigidBodies=True)
            phase3DumpDCD = dump.dcd(filename=self.outputDCD.replace('relaxed', 'phase3'), period=1, overwrite=True)
            phase3Step = integrate.mode_standard(dt = self.dtPhase3)
            phase3Flex = integrate.nvt(group=self.sideChainsGroup, T=self.T, tau=self.tau)
            phase3Rig = integrate.nvt_rigid(group=self.thioGroup, T=self.T, tau=self.tau)
            run(self.phase3RunLength)
            phase3DumpXML = dump.xml(filename=self.outputXML.replace('relaxed', 'phase3'), position = True, image = True, type = True, mass = True, diameter = True, body = True, charge = True, bond = True, angle = True, dihedral = True, improper = True)
            phase3Flex.disable()
            phase3Rig.disable()
            phase3DumpDCD.disable()
            del self.system, self.thioGroup, self.sideChainsGroup, self.energyLog, self.pair, self.b, self.a, self.d, self.i, phase3DumpDCD, phase3Step, phase3DumpXML, phase3Flex, phase3Rig
            init.reset()
        else:
            print "Phase 3 already completed for this morphology...skipping"
        # initDump.disable()
        # debugDump = dump.dcd(filename=self.outputDCD, period=1, overwrite=False)


        if self.runPhase4 == True:
            self.initialiseRun(self.outputXML.replace('relaxed', 'phase3'), pairType='lj', rigidBodies=True)
            phase4DumpDCD = dump.dcd(filename=self.outputDCD.replace('relaxed', 'phase4'), period=1, overwrite=True)
            phase4Step = integrate.mode_standard(dt=self.dtPhase4)
            phase4Flex = integrate.nvt(group=self.sideChainsGroup, T=self.T, tau=self.tau)
            phase4Rig = integrate.nvt_rigid(group=self.thioGroup, T=self.T, tau=self.tau)
            run(self.phase4RunLength)
            # self.initialPotentialEnergies = []
            # self.initialKineticEnergies = []
            # self.loadFromSnapshot = False
            # self.lowestKE = 9e999
            # self.KEIncreased = 0
            # self.firstKEValue = True
            # Modeler_hoomd wipes out the velocity data, and so we start with 0 kinetic energy. Skip this step (firstKEValue == True). Then, take the snapshot with the lowest KE.
            #checkKEs = analyze.callback(callback = self.checkKE, period=self.dumpPeriod)
            # try:
            #     run_upto(self.phase4RunLength)
            #     #run(self.maximumInitialRunLength)
            # except ExitHoomd as exitMessage:
            #     print exitMessage
            # if self.loadFromSnapshot == True:
            #     print "Loading from snapshot..."
            #     self.system.restore_snapshot(self.snapshotToLoad)
            phase4DumpXML = dump.xml(filename=self.outputXML.replace('relaxed', 'phase4'), position = True, image = True, type = True, mass = True, diameter = True, body = True, charge = True, bond = True, angle = True, dihedral = True, improper = True)
            phase4Flex.disable()
            phase4Rig.disable()
            phase4DumpDCD.disable()
            #checkKEs.disable()
            del self.system, self.thioGroup, self.sideChainsGroup, self.energyLog, self.pair, self.b, self.a, self.d, self.i, phase4DumpDCD, phase4Step, phase4Rig, phase4Flex, phase4DumpXML#, checkKEs, self.initialKineticEnergies, self.initialPotentialEnergies, self.loadFromSnapshot, exitMessage, self.snapshotToLoad
            # print dir()
            # print dir(self)
            init.reset()
        else:
            print "Phase 4 already completed for this morphology. Skipping..."


        if self.runPhase5 == True:
            self.initialiseRun(self.outputXML.replace('relaxed', 'phase4'), pairType='lj', rigidBodies = True)
            phase5DumpDCD = dump.dcd(filename=self.outputDCD.replace('relaxed', 'phase5'), period=1e3, overwrite=True)
            phase5Step = integrate.mode_standard(dt = self.dtPhase5)
            phase5Flex = integrate.nvt(group=self.sideChainsGroup, T=self.T, tau=self.tau)
            phase5Rig = integrate.nvt_rigid(group=self.thioGroup, T=self.T, tau=self.tau)
            run(self.phase5RunLength)
            phase5DumpXML = dump.xml(filename=self.outputXML.replace('relaxed', 'phase5'), position = True, image = True, type = True, mass = True, diameter = True, body = True, charge = True, bond = True, angle = True, dihedral = True, improper = True)
            phase5Flex.disable()
            phase5Rig.disable()
            phase5DumpDCD.disable()
            del self.system, self.thioGroup, self.sideChainsGroup, self.energyLog, self.pair, self.b, self.a, self.d, self.i, phase5DumpDCD, phase5Step, phase5Flex, phase5Rig, phase5DumpXML
            init.reset()
        else:
            print "Phase 5 already completed for this morphology...skipping"
        # initDump.disable()
        # debugDump = dump.dcd(filename=self.outputDCD, period=1, overwrite=False)




        ##### TESTING: FOR NOW JUST RUN THIS FOR AS LONG AS POSSIBLE TO SEE HOW THE ENERGY EVOLVES
        # return self.outputXML
        print os.listdir(self.saveDirectory)
        if (self.runPhase6 == True) or (self.continuePhase6 == True):
            # Then lock the sidechains in place and run the thiophenes for longer to make sure they equilibrate properly
            if self.continuePhase6 == False:
                self.initialiseRun(self.outputXML.replace('relaxed', 'phase5'), pairType = 'lj', rigidBodies = True)
            else:
                print "Continuing from previous run..."
                self.initialiseRun(self.continueFile)
                # underscoreList = helperFunctions.findIndex(self.continueFile, '_')
                # timestepsCompleted = int(self.continueFile[underscoreList[0]+1:underscoreList[1]])
                # self.mainRunLength -= timestepsCompleted
            phase6DumpDCD = dump.dcd(filename=self.outputDCD, period=self.mainTrajDumpPeriod, overwrite=True)
            phase6Step = integrate.mode_standard(dt=self.dtPhase6)
            phase6Flex = integrate.nvt(group=self.sideChainsGroup, T=self.T, tau=self.tau)
            phase6Rig = integrate.nvt_rigid(group=self.thioGroup, T=self.T, tau=self.tau)
            # self.mainPotentialEnergies = []
            # self.mainKineticEnergies = []
            # self.mainTotalEnergies = []
            # self.standardDeviation = []
            # self.maxStandardDeviation = 0
            # self.consecutiveDumpPeriodsUnderTarget = 0
            # checkTotalEs = analyze.callback(callback = self.getEnergies, period=self.dumpPeriod)
            resetXML = dump.xml(filename=self.outputXML.replace('relaxed_', 'temp_'), position = True, image = True, type = True, mass = True, diameter = True, body = True, charge = True, bond = True, angle = True, dihedral = True, improper = True, restart=True, period=self.phase6RunLength/10)
            try:
                run(self.phase6RunLength)
            except ExitHoomd as exitMessage:
                print exitMessage
            phase6DumpXML = dump.xml(filename=self.outputXML, position = True, image = True, type = True, mass = True, diameter = True, body = True, charge = True, bond = True, angle = True, dihedral = True, improper = True)
            phase6Flex.disable()
            phase6Rig.disable()
            phase6DumpDCD.disable()
            # checkTotalEs.disable()
            del self.system, self.thioGroup, self.sideChainsGroup, self.energyLog, self.pair, self.b, self.a, self.d, self.i, phase6DumpDCD, phase6Step, phase6Flex, phase6Rig, phase6DumpXML
            init.reset()
            # Run is complete so delete any temporary files
            saveFiles = os.listdir(self.saveDirectory)
            for fileName in saveFiles:
                if ('temp' in fileName) and (self.morphologyName in fileName):
                    print "Deleting temporary file:", self.saveDirectory+'/'+fileName
                    os.unlink(self.saveDirectory+'/'+fileName)
            # Plot PE
        else:
            print "Phase 6 already completed for this morphology. Skipping..."
        return self.outputXML

    
    def getEnergies(self, timestepNumber):
        currentPE = self.energyLog.query('potential_energy')
        currentKE = self.energyLog.query('kinetic_energy')
        currentPair = self.energyLog.query('pair_lj_energy')
        currentBond = self.energyLog.query('bond_harmonic_energy')
        currentAngle = self.energyLog.query('angle_harmonic_energy')
        currentDihedral = self.energyLog.query('dihedral_table_energy')
        self.mainPotentialEnergies.append(currentPE)
        self.mainKineticEnergies.append(currentKE)
        self.mainTotalEnergies.append(currentPE+currentKE+currentPair+currentBond+currentAngle+currentDihedral)
        currentStd = np.std(self.mainTotalEnergies[-10:])
        # print "Maximum STD =", self.maxStandardDeviation, "Current STD =", currentStd, "Target STD =", 0.05*self.maxStandardDeviation
        if currentStd > self.maxStandardDeviation:
            self.maxStandardDeviation = currentStd
            stdIncreasing = True
            self.consecutiveDumpPeriodsUnderTarget = 0
        else:
            stdIncreasing = False
        self.standardDeviation.append(currentStd)
        if (stdIncreasing == False) and (len(self.standardDeviation) >= 10):
            if currentStd <= 20:#0.05*self.maxStandardDeviation:
                self.consecutiveDumpPeriodsUnderTarget += 1
            else:
                self.consecutiveDumpPeriodsUnderTarget = 0
            if self.consecutiveDumpPeriodsUnderTarget == 10:
                raise ExitHoomd("Standard Deviation Condition Met", self.morphologyName)
        return 0

    def checkConvergence(self, timestepNumber):
        print "Flex Convergence =", self.phase2Flex.has_converged()
        #print "Rig Convergence =", self.phase2Rig.has_converged()
        if self.phase2Flex.has_converged():# and self.phase2Rig.has_converged():
            raise ExitHoomd("Position Minimisation has converged.")
            
    
        
    def checkKE(self, timestepNumber):
        currentPE = self.energyLog.query('potential_energy')
        currentKE = self.energyLog.query('kinetic_energy')
        if self.firstKEValue == False:
            if currentKE >= self.lowestKE:
                if self.KEIncreased == 5:
                    # Found the lowest KE point for at least 5 timesteps
                    del self.firstKEValue, self.lowestKE, self.KEIncreased
                    raise ExitHoomd("Lowest Energy Condition Met", self.morphologyName)
                self.KEIncreased += 1
            else:
                # Maybe at the lowest KE point so store snapshot
                self.KEIncreased = 0
                self.loadFromSnapshot = True
                self.snapshotToLoad = self.system.take_snapshot(all=True)
                self.lowestKE = currentKE
        else:
            self.firstKEValue = False
        self.initialPotentialEnergies.append(currentPE)
        self.initialKineticEnergies.append(currentKE)        
        del currentPE, currentKE
        return 0
        
    def fixCOM(self, atoms):
        # 'atoms' comes in as a list:
        # [[[functionGroup1AtomID1, functionalGroup1AtomID2,...], [targetCOM1]],...]
        for functionalGroup in atoms:
            print functionalGroup
            targetCOM = np.array(functionalGroup[1])
            currentCOM = self.calcCOM(functionalGroup[0])
            if (helperFunctions.findMagnitude(currentCOM-targetCOM) >= 3): #CURRENTLY IN ANGSTROEMS
                translation = targetCOM - currentCOM
                self.moveAtoms(functionalGroup[0], translation)

    def moveAtoms(self, atomIDs, translation):
        for atomID in atomIDs:
            currentPos = np.array(self.system.particles[atomID].position)
            newPos = list(currentPos + translation)
            self.system.particles[atomID].position = newPos

            
    def calcCOM(self, atomIDs):
        massWeightedX = 0.
        massWeightedY = 0.
        massWeightedZ = 0.
        totalMass = 0.
        for atomID in atomIDs:
            massWeightedX += self.system.particles[atomID].position[0]*self.system.particles[atomID].mass
            massWeightedY += self.system.particles[atomID].position[1]*self.system.particles[atomID].mass
            massWeightedZ += self.system.particles[atomID].position[2]*self.system.particles[atomID].mass
            totalMass += self.system.particles[atomID].mass
        #     print self.system.particles[atomID]
        #     print self.system.particles[atomID].position[0]
        # raise SystemError('STOP')
        return np.array([massWeightedX/float(totalMass), massWeightedY/float(totalMass), massWeightedZ/float(totalMass)])


    def plotEnergy(self, energies, fileName):
        timesteps = []
        for i in range(len(energies)):
            timesteps.append(i*self.dumpPeriod)
        plt.clf()
        plt.plot(timesteps, energies, 'ro')
        plt.ylabel('Energy')
        plt.xlabel('Timestep')
        plt.savefig(fileName)
    

def checkSaveDirectory(morphologyName, saveDirectory):
    saveDirectoryFiles = os.listdir(saveDirectory)
    runPhase1 = True
    runPhase2 = True
    runPhase3 = True
    runPhase4 = True
    runPhase5 = True
    runPhase6 = True
    continuePhase6 = False
    continueFile = None
    for fileName in saveDirectoryFiles:
        if morphologyName in fileName:
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
            elif ('temp' in fileName) and ('xml' in fileName):
                runPhase6 = False
                continuePhase6 = True
                continueFile = saveDirectory+fileName
    return [runPhase1, runPhase2, runPhase3, runPhase4, runPhase5, runPhase6, continuePhase6, continueFile]


def execute(morphologyFile, AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize):
    morphologyName = morphologyFile[helperFunctions.findIndex(morphologyFile,'/')[-1]+1:]
    outputDir = './outputFiles'
    morphologyList = os.listdir(outputDir)
    scaleFound = False
    for allMorphologies in morphologyList:
        if morphologyName in allMorphologies:
            outputDir += '/'+morphologyName
            break
    saveDir = outputDir+'/morphology'
    fileList = os.listdir(saveDir)
    continueData = checkSaveDirectory(morphologyName, saveDir)
    if np.sum(continueData[:-2]) == 0:
        exit()
    for fileName in fileList:
        if ('scaled' in fileName):
            scaleFound = True
            break
    # eScale = 1.
    # sScale = 1.
    eScale = 1./0.25
    sScale = 1./3.55
    slashList = helperFunctions.findIndex(AAfileName, '/')
    adjustedInputFileName = AAfileName[:slashList[-1]+1]+'scaled_'+str(1/sScale)+'_'+AAfileName[slashList[-1]+1:]

    if scaleFound == False:
        # Scale the positions relative to the origin
        # NOT ABLE TO USE MODELER HOOMD BECAUSE IT DOESN'T PRESERVE IMAGE LOCATIONS
        # WHICH IS VITAL WHEN RUNNING THE MORPHOLOGY IN ONE GO
        print "Scaling morphology..."
        if (sScale != 1.):
            AAMorphologyDict = helperFunctions.scale(AAMorphologyDict, sScale)
        print "Writing scaled XML..."
        helperFunctions.writeMorphologyXML(AAMorphologyDict, adjustedInputFileName)
    relaxedXML = hoomdRun(adjustedInputFileName, CGMoleculeDict, CGtoAAIDs, eScale, sScale, continueData, AAMorphologyDict).optimiseStructure()
    return morphologyFile, AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize
    

def loadPickle(morphologyFile):
    morphologyName = morphologyFile[helperFunctions.findIndex(morphologyFile,'/')[-1]+1:]
    outputDir = './outputFiles'
    morphologyList = os.listdir(outputDir)
    pickleFound = False
    for allMorphologies in morphologyList:
        if morphologyName in allMorphologies:
            outputDir += '/'+morphologyName
            break
    saveDir = outputDir+'/morphology'
    fileList = os.listdir(saveDir)
    for fileName in fileList:
        if fileName == morphologyName+'.pickle':
            pickleLoc = outputDir+'/morphology/'+fileName
            pickleFound = True
            break
    if pickleFound == False:
        print "Pickle file not found. Please run morphCT.py again to create the required HOOMD inputs."
        exit()
    print "Pickle found at", str(pickleLoc)+"."
    print "Loading data..."
    with open(pickleLoc, 'r') as pickleFile:
        (AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize) = pickle.load(pickleFile)
    continueData = checkSaveDirectory(morphologyName, saveDir)
    if np.sum(continueData[:-2]) != 0:
        morphologyFile, AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize = execute(morphologyFile, AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize)
    return morphologyFile, AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize
