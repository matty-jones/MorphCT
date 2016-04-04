from hoomd_script import *
import numpy as np
import modeler_hoomd as mh
import copy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
import helperFunctions
import sys

def multiHarmonicTorsion(theta, V0, V1, V2, V3, V4):
    #theta = np.deg2rad(theta)
    V = V0 + V1*np.cos(theta) + V2*((np.cos(theta))**2) + V3*((np.cos(theta))**3) + V4*((np.cos(theta))**4)
    F = V1*np.sin(theta) + 2*V2*np.cos(theta)*np.sin(theta) + 3*V3*((np.cos(theta))**2)*np.sin(theta) + 4*V4*((np.cos(theta))**3)*np.sin(theta)
    return (V, F)


class ExitHoomd(Exception):
    def __init__(self, string, moleculeName):
        self.string = string+" At Timestep = "+str(get_step())+" For Molecule = "+moleculeName
    def __str__(self):
        return self.string
    

class hoomdRun:
    def __init__(self, fileName, CGMoleculeDict, CG2AAIDs, eScale, sScale, continueData):
        slashLocs = helperFunctions.findIndex(fileName, '/')
        underscoreLocs = helperFunctions.findIndex(fileName, '_')
        self.morphologyName = fileName[underscoreLocs[-1]+1:-4]
        self.saveDirectory = fileName[:slashLocs[-1]]+'/'
        self.fileName = fileName
        self.CGMoleculeDict = CGMoleculeDict
        self.CG2AAIDs = CG2AAIDs
        self.eScale = eScale
        self.sScale = sScale
        self.dumpPeriod = 1e2
        self.mainTrajDumpPeriod = 1e3
        self.T = 1.0
        self.tau = 1.0
        self.dtPhase1 = 1e-5
        self.dtPhase2 = 1e-4
        self.dtPhase3 = 2e-4
        self.dtPhase4 = 1e-3
        self.dtPhase5 = 1e-3
        self.phase1RunLength = 1000
        self.phase2RunLength = 50000
        self.phase3RunLength = 50000
        self.phase4RunLength = 1e6 # This is a maximum because sim is curtailed
        self.phase5RunLength = 5e4
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
        self.continuePhase5 = continueData[5]
        self.continueFile = continueData[6]
        

    def initialiseRun(self, inputFilename):
        self.system = init.read_xml(filename=inputFilename)
        self.thioGroup = group.tag_list(name="thio", tags=self.thioGroupIDs)
        self.alk1Group = group.tag_list(name="alk1", tags=self.alk1GroupIDs)
        self.alk2Group = group.tag_list(name="alk2", tags=self.alk2GroupIDs)
        self.setForcefieldCoeffs(self.eScale, self.sScale)
        self.energyLog = analyze.log(filename = self.outputLOG, quantities=['potential_energy', 'kinetic_energy', 'pair_lj_energy', 'bond_harmonic_energy', 'angle_harmonic_energy', 'dihedral_table_energy'], period=self.dumpPeriod, overwrite = self.overwriteEnergies)
        self.overwriteEnergies = False


    def setForcefieldCoeffs(self, eScale, sScale):
        # Forcefield parameters obtained from:
        # Bhatta, R. S., Yimer, Y. Y, Perry, D. S., Tsige, M., "Improved Force Field for Molecular Modeling of Poly(3-Hexylthiophene)", 2013, J. Phys. Chem. B, DOI: 10.1021/jp404629a
        # Marcon, V., Raos, G., "Free Energies of Molecular Crystal Surfaces by Computer Simlations: Application to Tetrathiophene", 2006, J. Am. Chem. Soc., DOI: 10.1021/ja056548t
        # Marcon, V., Raos, G., "Molecular Modeling of Crystalline Oligothiophenes: Testing and Development of Improved Force Fields", 2004, J. Phys. Chem. B, DOI: 10.1021/jp047128d
        # Huang, D. M., Faller, R., Do, K., Moule, A. J., "Coarse-Grained Computer Simulations of Polymer/Fullerene Bulk Heterojunctions for Organic Photovoltaic Applications", 2010, J. Chem. Theory Comput., DOI: 10.1021/ct900496t
        # Jorgensen, W. L., Maxwell, D. S., Tirdao-Rives, J., "Development and Testing of the OPLS All-Atom Forcefield on Conformational Energetics and Properties of Organic Liquids", 1996, J. Am. Chem. Soc. DOI: 10.1021/ja9621760
        # Bhatta, R. S., Yimer, Y. Y., Tsige, M., Perry, D. S., "Conformations and Torsional Potentials of Poly(3-Hexylthiophene) Oligomers: Density Functional Calculations Up to the Dodecamer", 2012, Comput. & Theor. Chem., DOI: 10.1016/j.comptc.2012.06.026

        self.lj = pair.lj(r_cut=10*sScale)
        self.lj.set_params(mode='xplor')
        self.lj.pair_coeff.set('C1','C1',epsilon=0.070*eScale,sigma=3.550*sScale)
        self.lj.pair_coeff.set('C1','C2',epsilon=0.070*eScale,sigma=3.550*sScale)
        self.lj.pair_coeff.set('C1','C3',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C1','C4',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C1','C5',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C1','C6',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C1','C7',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C1','C8',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C1','C9',epsilon=0.070*eScale,sigma=3.550*sScale)
        self.lj.pair_coeff.set('C1','C10',epsilon=0.070*eScale,sigma=3.550*sScale)
        self.lj.pair_coeff.set('C1','H1',epsilon=0.046*eScale,sigma=2.979*sScale)
        self.lj.pair_coeff.set('C1','S1',epsilon=0.132*eScale,sigma=3.550*sScale)
        self.lj.pair_coeff.set('C2','C2',epsilon=0.250*eScale,sigma=3.550*sScale)
        self.lj.pair_coeff.set('C2','C3',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C2','C4',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C2','C5',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C2','C6',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C2','C7',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C2','C8',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C2','C9',epsilon=0.070*eScale,sigma=3.550*sScale)
        self.lj.pair_coeff.set('C2','C10',epsilon=0.070*eScale,sigma=3.550*sScale)
        self.lj.pair_coeff.set('C2','H1',epsilon=0.046*eScale,sigma=2.979*sScale)
        self.lj.pair_coeff.set('C2','S1',epsilon=0.132*eScale,sigma=3.550*sScale)
        self.lj.pair_coeff.set('C3','C3',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.lj.pair_coeff.set('C3','C4',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.lj.pair_coeff.set('C3','C5',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.lj.pair_coeff.set('C3','C6',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.lj.pair_coeff.set('C3','C7',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.lj.pair_coeff.set('C3','C8',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.lj.pair_coeff.set('C3','C9',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C3','C10',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C3','H1',epsilon=0.044*eScale,sigma=2.958*sScale)
        self.lj.pair_coeff.set('C3','S1',epsilon=0.128*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C4','C4',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.lj.pair_coeff.set('C4','C5',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.lj.pair_coeff.set('C4','C6',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.lj.pair_coeff.set('C4','C7',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.lj.pair_coeff.set('C4','C8',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.lj.pair_coeff.set('C4','C9',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C4','C10',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C4','H1',epsilon=0.044*eScale,sigma=2.958*sScale)
        self.lj.pair_coeff.set('C4','S1',epsilon=0.128*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C5','C5',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.lj.pair_coeff.set('C5','C6',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.lj.pair_coeff.set('C5','C7',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.lj.pair_coeff.set('C5','C8',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.lj.pair_coeff.set('C5','C9',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C5','C10',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C5','H1',epsilon=0.044*eScale,sigma=2.958*sScale)
        self.lj.pair_coeff.set('C5','S1',epsilon=0.128*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C6','C6',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.lj.pair_coeff.set('C6','C7',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.lj.pair_coeff.set('C6','C8',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.lj.pair_coeff.set('C6','C9',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C6','C10',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C6','H1',epsilon=0.044*eScale,sigma=2.958*sScale)
        self.lj.pair_coeff.set('C6','S1',epsilon=0.128*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C7','C7',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.lj.pair_coeff.set('C7','C8',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.lj.pair_coeff.set('C7','C9',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C7','C10',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C7','H1',epsilon=0.044*eScale,sigma=2.958*sScale)
        self.lj.pair_coeff.set('C7','S1',epsilon=0.128*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C8','C8',epsilon=0.066*eScale,sigma=3.500*sScale)
        self.lj.pair_coeff.set('C8','C9',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C8','C10',epsilon=0.068*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C8','H1',epsilon=0.044*eScale,sigma=2.958*sScale)
        self.lj.pair_coeff.set('C8','S1',epsilon=0.128*eScale,sigma=3.525*sScale)
        self.lj.pair_coeff.set('C9','C9',epsilon=0.070*eScale,sigma=3.550*sScale)
        self.lj.pair_coeff.set('C9','C10',epsilon=0.070*eScale,sigma=3.500*sScale)
        self.lj.pair_coeff.set('C9','H1',epsilon=0.046*eScale,sigma=2.979*sScale)
        self.lj.pair_coeff.set('C9','S1',epsilon=0.132*eScale,sigma=3.550*sScale)
        self.lj.pair_coeff.set('C10','C10',epsilon=0.070*eScale,sigma=3.550*sScale)
        self.lj.pair_coeff.set('C10','H1',epsilon=0.046*eScale,sigma=2.979*sScale)
        self.lj.pair_coeff.set('C10','S1',epsilon=0.132*eScale,sigma=3.550*sScale)
        self.lj.pair_coeff.set('H1','H1',epsilon=0.030*eScale,sigma=2.500*sScale)
        self.lj.pair_coeff.set('H1','S1',epsilon=0.087*eScale,sigma=2.979*sScale)
        self.lj.pair_coeff.set('S1','S1',epsilon=0.250*eScale,sigma=3.550*sScale)

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

        self.a = angle.harmonic()
        # Ring Bond Angles [k] = kcal mol^{-1} rad^{-2} * epsilon, [t] = rad
        self.a.set_coeff('C10-C9-C2',k=129940*eScale,t0=1.97784)
        self.a.set_coeff('C10-C9-H1',k=115762*eScale,t0=2.14639)
        self.a.set_coeff('C10-S1-C1',k=283503*eScale,t0=1.61921)
        self.a.set_coeff('C9-C10-S1',k=283503*eScale,t0=1.92496)
        self.a.set_coeff('C1-C2-C9',k=129940*eScale,t0=1.97784)
        self.a.set_coeff('C2-C9-H1',k=115762*eScale,t0=2.15897)
        self.a.set_coeff('C2-C1-S1',k=283503*eScale,t0=1.92496)
        # Alkyl Bond Angles
        self.a.set_coeff('C3-C2-C9',k=546735*eScale,t0=2.15335)
        self.a.set_coeff('C2-C3-C4',k=394396*eScale,t0=2.01481)
        self.a.set_coeff('C2-C3-H1',k=243125*eScale,t0=1.90571)
        self.a.set_coeff('C1-C2-C3',k=545996*eScale,t0=2.17388)
        self.a.set_coeff('C3-C4-C5',k=191552*eScale,t0=1.96699)
        self.a.set_coeff('C3-C4-H1',k=123105*eScale,t0=1.93208)
        self.a.set_coeff('C4-C3-H1',k=123105*eScale,t0=1.93208)
        self.a.set_coeff('C4-C5-C6',k=191552*eScale,t0=1.96699)
        self.a.set_coeff('C4-C5-H1',k=123105*eScale,t0=1.93208)
        self.a.set_coeff('C5-C4-H1',k=123105*eScale,t0=1.93208)
        self.a.set_coeff('C5-C6-C7',k=191552*eScale,t0=1.96699)
        self.a.set_coeff('C5-C6-H1',k=123105*eScale,t0=1.93208)
        self.a.set_coeff('C6-C5-H1',k=123105*eScale,t0=1.93208)
        self.a.set_coeff('C6-C7-C8',k=191552*eScale,t0=1.96699)
        self.a.set_coeff('C6-C7-H1',k=123105*eScale,t0=1.93208)
        self.a.set_coeff('C7-C6-H1',k=123105*eScale,t0=1.93208)
        self.a.set_coeff('C7-C8-H1',k=123105*eScale,t0=1.93208)
        self.a.set_coeff('C8-C7-H1',k=123105*eScale,t0=1.93208)
        self.a.set_coeff('H1-C8-H1',k=108333*eScale,t0=1.88146)
        self.a.set_coeff('H1-C7-H1',k=108333*eScale,t0=1.88146)
        self.a.set_coeff('H1-C6-H1',k=108333*eScale,t0=1.88146)
        self.a.set_coeff('H1-C5-H1',k=108333*eScale,t0=1.88146)
        self.a.set_coeff('H1-C4-H1',k=108333*eScale,t0=1.88146)
        self.a.set_coeff('H1-C3-H1',k=108333*eScale,t0=1.88146)
        # Inter-monomer Bond Angles
        self.a.set_coeff('C1-C10-C9',k=179550*eScale,t0=2.27137)
        self.a.set_coeff('C1-C10-S1',k=137024*eScale,t0=2.08687)
        self.a.set_coeff('S1-C1-C10',k=137024*eScale,t0=2.08687)
        self.a.set_coeff('C2-C1-C10',k=179550*eScale,t0=2.27137)

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
                if self.CG2AAIDs[moleculeNo][CGAtomID][0] == 'thio':
                    thioAtoms.append(self.CG2AAIDs[moleculeNo][CGAtomID][1:])
                    self.thioGroupIDs += self.CG2AAIDs[moleculeNo][CGAtomID][1]
                elif self.CG2AAIDs[moleculeNo][CGAtomID][0] == 'alk1':
                    alk1Atoms.append(self.CG2AAIDs[moleculeNo][CGAtomID][1:])
                    self.alk1GroupIDs += self.CG2AAIDs[moleculeNo][CGAtomID][1]
                elif self.CG2AAIDs[moleculeNo][CGAtomID][0] == 'alk2':
                    alk2Atoms.append(self.CG2AAIDs[moleculeNo][CGAtomID][1:])
                    self.alk2GroupIDs += self.CG2AAIDs[moleculeNo][CGAtomID][1]
        # Build the groups:
        if self.runPhase1 == True:
            self.initialiseRun(self.fileName)
            phase1DumpDCD = dump.dcd(filename=self.outputDCD.replace('relaxed', 'phase1'), period=10, overwrite=True)
            phase1Step = integrate.mode_standard(dt = self.dtPhase1)
            # phase1 = integrate.brownian(group=group.all(), seed=3, dscale=1e11, T=self.T)
            phase1 = integrate.nve(group=group.all(), limit=0.01)
            run(self.phase1RunLength)
            phase1DumpXML = dump.xml(filename=self.outputXML.replace('relaxed', 'phase1'), all=True)
            phase1.disable()
            phase1DumpDCD.disable()
            del self.system, self.thioGroup, self.alk1Group, self.alk2Group, self.energyLog, self.lj, self.b, self.a, self.d, self.i, phase1DumpDCD, phase1Step, phase1, phase1DumpXML
            init.reset()
        else:
            print "Phase 1 already completed for this morphology...skipping"

            
        if self.runPhase2 == True:
            self.initialiseRun(self.outputXML.replace('relaxed', 'phase1'))
            phase2DumpDCD = dump.dcd(filename=self.outputDCD.replace('relaxed', 'phase2'), period=100, overwrite=True)
            phase2Step = integrate.mode_standard(dt = self.dtPhase2)
            phase2 = integrate.nvt(group=group.all(), T=self.T, tau=self.tau)
            run(self.phase2RunLength)
            phase2DumpXML = dump.xml(filename=self.outputXML.replace('relaxed', 'phase2'), all=True)
            phase2.disable()
            phase2DumpDCD.disable()
            del self.system, self.thioGroup, self.alk1Group, self.alk2Group, self.energyLog, self.lj, self.b, self.a, self.d, self.i, phase2DumpDCD, phase2Step, phase2, phase2DumpXML
            init.reset()
        else:
            print "Phase 2 already completed for this morphology...skipping"
        # initDump.disable()
        # debugDump = dump.dcd(filename=self.outputDCD, period=1, overwrite=False)

        if self.runPhase3 == True:
            self.initialiseRun(self.outputXML.replace('relaxed', 'phase2'))
            phase3DumpDCD = dump.dcd(filename=self.outputDCD.replace('relaxed', 'phase3'), period=100, overwrite=True)
            phase3Step = integrate.mode_standard(dt = self.dtPhase3)
            phase3 = integrate.nvt(group=group.all(), T=self.T, tau=self.tau)
            run(self.phase3RunLength)
            phase3DumpXML = dump.xml(filename=self.outputXML.replace('relaxed', 'phase3'), all=True)
            phase3.disable()
            phase3DumpDCD.disable()
            del self.system, self.thioGroup, self.alk1Group, self.alk2Group, self.energyLog, self.lj, self.b, self.a, self.d, self.i, phase3DumpDCD, phase3Step, phase3, phase3DumpXML
            init.reset()
        else:
            print "Phase 3 already completed for this morphology...skipping"
        # initDump.disable()
        # debugDump = dump.dcd(filename=self.outputDCD, period=1, overwrite=False)


        if self.runPhase4 == True:
            self.initialiseRun(self.outputXML.replace('relaxed', 'phase3'))
            phase4DumpDCD = dump.dcd(filename=self.outputDCD.replace('relaxed', 'phase4'), period=self.dumpPeriod, overwrite=True)
            phase4Step = integrate.mode_standard(dt=self.dtPhase4)
            phase4 = integrate.nvt(group=group.all(), T=self.T, tau=self.tau)
            # timesteps = 6
            # currentTimestep = 1
            # while True:
            #     if currentTimestep == timesteps:
            #         break
            #     # snapshot = moleculeData.take_snapshot()
            #     # LEAVE COM FIX FOR NOW, INSTEAD LOCK SIDECHAINS AFTER SHORT TRAJ RUN
            #     # self.fixCOM(thioAtoms)
            #     # self.fixCOM(alk1Atoms)
            #     # self.fixCOM(alk2Atoms)
            #     run(1e3)
            #     currentTimestep += 1
            self.initialPotentialEnergies = []
            self.initialKineticEnergies = []
            self.loadFromSnapshot = False
            self.lowestKE = 9e999
            self.KEIncreased = 0
            self.firstKEValue = True
            # Modeler_hoomd wipes out the velocity data, and so we start with 0 kinetic energy. Skip this step (firstKEValue == True). Then, take the snapshot with the lowest KE.
            checkKEs = analyze.callback(callback = self.checkKE, period=self.dumpPeriod)
            try:
                run_upto(self.phase4RunLength)
                #run(self.maximumInitialRunLength)
            except ExitHoomd as exitMessage:
                print exitMessage
            if self.loadFromSnapshot == True:
                print "Loading from snapshot..."
                self.system.restore_snapshot(self.snapshotToLoad)
            phase4DumpXML = dump.xml(filename=self.outputXML.replace('relaxed', 'phase3'), all=True)
            phase4.disable()
            phase4DumpDCD.disable()
            checkKEs.disable()
            del self.system, self.snapshotToLoad, self.thioGroup, self.alk1Group, self.alk2Group, self.energyLog, self.lj, self.b, self.a, self.d, self.i, phase4DumpDCD, phase4Step, phase4, phase4DumpXML, checkKEs, self.initialKineticEnergies, self.initialPotentialEnergies, self.loadFromSnapshot, exitMessage
            # print dir()
            # print dir(self)
            init.reset()
        else:
            print "Phase 4 already completed for this morphology. Skipping..."

            
        ##### TESTING: FOR NOW JUST RUN THIS FOR AS LONG AS POSSIBLE TO SEE HOW THE ENERGY EVOLVES
        # return self.outputXML
        print os.listdir(self.saveDirectory)
        if (self.runPhase5 == True) or (self.continuePhase5 == True):
            # Then lock the sidechains in place and run the thiophenes for longer to make sure they equilibrate properly
            if self.continuePhase5 == False:
                self.initialiseRun(self.outputXML.replace('relaxed', 'phase4'))
            else:
                self.initialiseRun(self.continueFile)
                # underscoreList = helperFunctions.findIndex(self.continueFile, '_')
                # timestepsCompleted = int(self.continueFile[underscoreList[0]+1:underscoreList[1]])
                # self.mainRunLength -= timestepsCompleted
            phase5DumpDCD = dump.dcd(filename=self.outputDCD, period=self.mainTrajDumpPeriod, overwrite=True)
            phase5Step = integrate.mode_standard(dt=self.dtPhase5)
            phase5 = integrate.nvt(group=self.thioGroup, T=self.T, tau=self.tau)
            # self.mainPotentialEnergies = []
            # self.mainKineticEnergies = []
            # self.mainTotalEnergies = []
            # self.standardDeviation = []
            # self.maxStandardDeviation = 0
            # self.consecutiveDumpPeriodsUnderTarget = 0
            # checkTotalEs = analyze.callback(callback = self.getEnergies, period=self.dumpPeriod)
            resetXML = dump.xml(filename=self.outputXML.replace('relaxed_', 'temp_'), all=True, restart=True, period=self.mainRunLength/10)
            try:
                run_upto(self.phase5RunLength)
            except ExitHoomd as exitMessage:
                print exitMessage
            phase5DumpXML = dump.xml(filename=self.outputXML, all=True)
            phase5.disable()
            phase5DumpDCD.disable()
            # checkTotalEs.disable()
            del self.system, self.thioGroup, self.alk1Group, self.alk2Group, self.energyLog, self.lj, self.b, self.a, self.d, self.i, phase5DumpDCD, phase5Step, phase5, phase5DumpXML
            init.reset()
            # Run is complete so delete any temporary files
            saveFiles = os.listdir(self.saveDirectory)
            for fileName in saveFiles:
                if ('temp' in fileName) and (self.morphologyName in fileName):
                    print "Deleting temporary file:", self.saveDirectory+'/'+fileName
                    os.unlink(self.saveDirectory+'/'+fileName)
            # Plot PE
        else:
            print "Phase 5 already completed for this morphology. Skipping..."
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
    continuePhase5 = False
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
            elif ('temp' in fileName) and ('xml' in fileName):
                runPhase5 = False
                continuePhase5 = True
                continueFile = saveDirectory+fileName
    return [runPhase1, runPhase2, runPhase3, runPhase4, runPhase5, continuePhase5, continueFile]
        
#def execute(moleculeDir):
if __name__ == '__main__':
    morphologyFile = sys.argv[1]
    morphologyName = morphologyFile[helperFunctions.findIndex(morphologyFile,'/')[-1]+1:]
    outputDir = './outputFiles'
    morphologyList = os.listdir(outputDir)
    pickleFound = False
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
        if fileName == morphologyName+'.pickle':
            pickleLoc = outputDir+'/morphology/'+fileName
            pickleFound = True
    if pickleFound == False:
        print "Pickle file not found. Please run morphCT.py again to create the required HOOMD inputs."
        exit()
    print "Pickle found at", str(pickleLoc)+"."
    print "Loading data..."
    with open(pickleLoc, 'r') as pickleFile:
        (AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, boxSize) = pickle.load(pickleFile)
    eScale = 1.
    sScale = 1.
    #eScale = 1./0.25
    #sScale = 1./3.55
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

    # # Rescale morphologies using modeler_hoomd
    # buildingBlock = mh.builder.building_block(AAfileName, model_name = '3HTmid', model_file = os.getcwd()+'/templates/model.xml')
    # print "Scaling morphology..."
    # if (sScale != 1.):
    #     buildingBlock.scale(sScale)
    # print "Writing scaled XML..."
    # buildingBlock.write(adjustedInputFileName)
    # print "Updating simulation dimensions..."
    # helperFunctions.updateXMLBoxLength(adjustedInputFileName, boxSize)

    relaxedXML = hoomdRun(adjustedInputFileName, CGMoleculeDict, CGtoAAIDs, eScale, sScale, continueData).optimiseStructure()

    # raise SystemError('HALT')
    # # Now make the .xyz for Lan
    # # First rescale
    # xyzName = relaxedXML[:-3]+'xyz'

    # print "Run Complete."

    # print "Rescaling morphology..."
    # buildingBlock = mh.builder.building_block(relaxedXML, model_name = '3HT_mid', model_file = os.getcwd()+'/templates/model.xml')
    # if (sScale != 1.):
    #     buildingBlock.scale(1/sScale)
    # print "Writing rescaled temporary XML..."
    # buildingBlock.write(os.getcwd()+'/temporary_molecule.xml')
    # print "Converting to .xyz for DFT..."
    # helperFunctions.XMLToXYZ(os.getcwd()+'/temporary_molecule.xml', xyzName)
    # print "Tidying up..."
    # os.unlink(os.getcwd()+'/temporary_molecule.xml')
    # init.reset()
