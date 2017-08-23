import runMorphCT

# ---=== Directory and File Structure ===---
inputDir = '/Users/mattyjones/GoogleDrive/Boise/Code/MorphCT/inputCGMorphs'
outputDir = '/Users/mattyjones/GoogleDrive/Boise/Code/MorphCT/outputFiles'


# ---=== Input Morphology Details ===---
morphology = 'peryleneBlendTest.xml'
inputSigma = 1.0
overwriteCurrentData = True

# ---=== Execution Modules ===---

executeFinegraining = True
executeMolecularDynamics = True
executeExtractMolecules = False
executeObtainChromophores = False
executeZINDO = False
executeCalculateTransferIntegrals = False
executeCalculateMobility = False

# ---=== Fine Graining Parameters ===---

CGToTemplateDirs = {\
'A':'/Users/mattyjones/GoogleDrive/Boise/Code/MorphCT/templates',\
'B':'/Users/mattyjones/GoogleDrive/Boise/Code/MorphCT/templates',\
}
CGToTemplateFiles = {\
'A':'perylene.xml',\
'B':'perylothiophene.xml',\
}
CGToTemplateForceFields = {\
'A':'FFPerylene.xml',\
'B':'FFPerylothiophene.xml',\
}
CGToTemplateAAIDs = {\
'A':[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31],\
'B':[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30],\
}
CGToTemplateBonds = {\
}
rigidBodySites = {\
#'A':[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31],\
#'B':[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30],\
}
additionalConstraints = [\
]
moleculeTerminatingUnits = {\
}
moleculeTerminatingBonds = [\
]
moleculeTerminatingConnections = [\
]

# ---=== Forcefield Parameters ===---
pairRCut = 10
pairDPDGammaVal = 0.0
## --== Lennard-Jones Pair ==--
#ljCoeffs = [\
## Perylene
#['C30-C30', 0.0, 0.0],
#['C30-C31', 0.0, 0.0],
#['C30-H30', 0.0, 0.0],
## Perylothiophene
#['C40-C40', 0.0, 0.0],
#['C40-C41', 0.0, 0.0],
#['C30-H30', 0.0, 0.0],
#]
## --== Disipative Particle Dynamics Pair ==--
#dpdCoeffs = [\
#]
## --== Bond ==--
#bondCoeffs = [\
#]
## --== Angle ==--
#angleCoeffs = [\
#]
## --== Dihedral ==--
#dihedralCoeffs = [\
#]
## --== Improper ==--
#improperCoeffs = [\
#]

# ---=== Molecular Dynamics Phase Parameters ===---
numberOfPhases = 8
temperatures = [1.0]
taus = [1.0]
pairTypes = ['none', 'dpd', 'lj', 'lj', 'lj', 'lj', 'lj', 'lj']
bondTypes = ['harmonic']
angleTypes = ['harmonic']
dihedralTypes = ['opls']
integrationTargets = ['all']
timesteps = [1E-3, 1E-3, 1E-9, 5E-9, 1E-8, 1E-7, 1E-6, 1E-5]
durations = [1E5, 1E4, 1E3, 1E3, 1E3, 1E4, 1E5, 1E5]
terminationConditions = ['KEmin', 'maxt', 'maxt', 'maxt', 'maxt', 'maxt', 'maxt', 'maxt']
groupAnchorings = ['all', 'all', 'all', 'all', 'all', 'all', 'all', 'none']

# ---=== Chromophore Parameters ===---

CGSiteSpecies = {\
'A':'Donor',\
'B':'None',\
'C':'None',\
'D':'Acceptor',\
}
maximumHopDistance = 10.0
removeORCAInputs = False
removeORCAOutputs = False
chromophoreLength = 3

# ---=== Chromophore Energy Scaling Parameters ===---

literatureHOMO = -5.0
literatureLUMO = None
targetDoSSTDHOMO = 0.1
targetDoSSTDLUMO = None


# ---=== Kinetic Monte Carlo Parameters ===---

systemTemperature = 290
numberOfCarriersPerSimulationTime = 10
hopLimit = 0
simulationTimes = [1.00e-6]
recordCarrierHistory = True
reorganisationEnergy = 0.3063
combineKMCResults = False


# ---=== Begin run ===---
parameterFile = __file__

if __name__ == "__main__":
    parameterNames = [i for i in dir() if (not i.startswith('__')) and (i not in ['runMorphCT', 'os', 'sys'])]
    parameters = {}
    for name in parameterNames:
        parameters[name] = locals()[name]
    runMorphCT.simulation(**parameters) # Execute MorphCT using these simulation parameters
