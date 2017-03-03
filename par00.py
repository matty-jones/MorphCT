# ---=== Directory and File Structure ===---
inputDir = '/Users/mattyjones/GoogleDrive/Boise/Code/MorphCT/inputCGMorphs'
outputDir = '/Users/mattyjones/GoogleDrive/Boise/Code/MorphCT/outputFiles'


# ---=== Input Morphology Details ===---
morphology = 'p1-L15-f0.3-P0.1-T1.0-e0.1.xml'
inputSigma = 3.0
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
'C':'/Users/mattyjones/GoogleDrive/Boise/Code/MorphCT/templates',\
'D':'/Users/mattyjones/GoogleDrive/Boise/Code/MorphCT/templates',\
}
CGToTemplateFiles = {\
'A':'mid3HT.xml',\
'B':'mid3HT.xml',\
'C':'mid3HT.xml',\
'D':'PCBM.xml',\
}
CGToTemplateForceFields = {\
'A':'FFP3HT.xml',\
'B':'FFP3HT.xml',\
'C':'FFP3HT.xml',\
'D':'FFPCBM.xml',\
}
CGToTemplateAAIDs = {\
'A':[0, 1, 2, 3, 4, 24],\
'B':[5, 6, 7, 18, 19, 20, 21, 22, 23],\
'C':[8, 9, 10, 11, 12, 13, 14, 15, 16, 17],\
#'D':[],\
'D':[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87],\
}
CGToTemplateBonds = {\
'bondB':['C2-C3', 2, 5],\
'bondC':['C5-C6', 7, 8],\
}
rigidBodySites = {\
'A':[0, 1, 2, 3, 4],\
#'D':[],\
'D':[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59],\
}
additionalConstraints = [\
['C1-C10', 3, 25],\
['C1-C10-C9', 3, 25, 26],\
['C1-C10-S1', 3, 25, 29],\
['S1-C1-C10', 4, 3, 25],\
['C2-C1-C10', 2, 3, 25],\
['C1-C10-C9-C2', 3, 25, 26, 27],\
['C1-C10-S1-C1', 3, 25, 29, 28],\
['S1-C1-C10-S1', 4, 3, 25, 29],\
['C2-C1-C10-S1', 2, 3, 25, 29],\
]
moleculeTerminatingUnits = [\
['H1',0,0,0],\
]
moleculeTerminatingBonds = [\
]
moleculeTerminatingConnections = [\
['C10-H1', 0],\
['C1-H1', 3],\
]

# ---=== Forcefield Parameters ===---

pairRCut = 10.0
pairDPDGammaVal = 0.0

# ---=== Molecular Dynamics Phase Parameters ===---
numberOfPhases = 8
temperatures = [1.0]
taus = [1.0]
pairTypes = ['none', 'dpd', 'lj', 'lj', 'lj', 'lj', 'lj', 'lj']
bondTypes = ['harmonic']
angleTypes = ['harmonic']
dihedralTypes = ['table']
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
    import runMorphCT
    import sys

    sys.path.append('./code')

    import helperFunctions

    procIDs = helperFunctions.getCPUCores()
    parameterNames = [i for i in dir() if (not i.startswith('__')) and (i not in ['runMorphCT', 'helperFunctions', 'sys'])]
    parameters = {}
    for name in parameterNames:
        parameters[name] = locals()[name]
    runMorphCT.simulation(**parameters) # Execute MorphCT using these simulation parameters
