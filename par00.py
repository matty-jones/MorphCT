# ---=== Directory and File Structure ===---
inputDir = '/scratch/erjank_project/matty/blendTest/inputCGMorphs'
outputDir = '/scratch/erjank_project/matty/blendTest/outputFiles'

# ---=== Input Morphology Details ===---

morphology = 'P3HTPCBM.xml'
inputSigma = 3.0
overwriteCurrentData = True

# ---=== Execution Modules ===---

executeFinegraining = True
executeMolecularDynamics = True
executeExtractMolecules = False
executeObtainChromophores = True
executeZINDO = True
executeCalculateTransferIntegrals = True
executeCalculateMobility = True

# ---=== Fine Graining Parameters ===---

CGToTemplateDirs = {\
'A':'/scratch/erjank_project/matty/blendTest/templates',\
'B':'/scratch/erjank_project/matty/blendTest/templates',\
'C':'/scratch/erjank_project/matty/blendTest/templates',\
'D':'/scratch/erjank_project/matty/blendTest/templates',\
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
'D':[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87],\
}
CGToTemplateBonds = {\
'bondB':['C2-C3', 2, 5],\
'bondC':['C5-C6', 7, 8],\
}
rigidBodySites = {\
'A':[0, 1, 2, 3, 4],\
'D':[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87],\
}
additionalConstraints = [\
['A', 'C1-C10', 3, 25],\
['A', 'C1-C10-C9', 3, 25, 26],\
['A', 'C1-C10-S1', 3, 25, 29],\
['A', 'S1-C1-C10', 4, 3, 25],\
['A', 'C2-C1-C10', 2, 3, 25],\
['A', 'C1-C10-C9-C2', 3, 25, 26, 27],\
['A', 'C1-C10-S1-C1', 3, 25, 29, 28],\
['A', 'S1-C1-C10-S1', 4, 3, 25, 29],\
['A', 'C2-C1-C10-S1', 2, 3, 25, 29],\
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
timesteps = [1E-3]
durations = [1E3]
terminationConditions = ['maxt']
groupAnchorings = ['all']
DCDFileWrite = True
DCDFileDumpsteps = [0]

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
literatureLUMO = -4.5
targetDoSSTDHOMO = 0.1
targetDoSSTDLUMO = 0.1


# ---=== Kinetic Monte Carlo Parameters ===---

systemTemperature = 290
<<<<<<< Updated upstream
numberOfCarriersPerSimulationTime = 1000
=======
numberOfCarriersPerSimulationTime = 100
>>>>>>> Stashed changes
hopLimit = 0
simulationTimes = [1.00e-12, 1.00e-11, 1.00e-10, 1.00e-9, 1.00e-8]
recordCarrierHistory = True
reorganisationEnergyDonor = 0.3063
reorganisationEnergyAcceptor = 0.1600
combineKMCResults = True

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
