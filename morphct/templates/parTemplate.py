# ---=== Directory and File Structure ===---
inputDir = INPUTDIR
outputDir = OUTPUTDIR

# ---=== Input Morphology Details ===---

morphology = INPUTMORPHOLOGY
inputSigma = INPUTSIGMA
overwriteCurrentData = True

# ---=== Execution Modules ===---

executeFinegraining = True
executeMolecularDynamics = True
executeExtractMolecules = True
executeObtainChromophores = True
executeZINDO = True
executeCalculateTransferIntegrals = True
executeCalculateMobility = True

# ---=== Fine Graining Parameters ===---

CGToTemplateDirs = {\
CGTOTEMPLATEDIRS
}
CGToTemplateFiles = {\
CGTOTEMPLATEFILES
}
CGToTemplateForceFields = {\
CGTOTEMPLATEFORCEFIELDS
}
CGToTemplateAAIDs = {\
CGTOTEMPLATEAAIDS
}
CGToTemplateBonds = {\
CGTOTEMPLATEBONDS
}
rigidBodySites = {\
RIGIDBODYSITES
}
additionalConstraints = [\
ADDITIONALCONSTRAINTS
]
moleculeTerminatingUnits = [\
'H1': [0, 0, 0]
]
moleculeTerminatingBonds = [\
]
moleculeTerminatingConnections = [\
TERMINATINGCONNECTIONS
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
dihedralTypes = ['opls']
integrationTargets = ['all']
timesteps = [1E-3, 1E-3, 1E-10, 1E-9, 1E-8, 1E-7, 1E-6, 1E-5]
durations = [1E5, 1E4, 1E3, 1E3, 1E3, 1E4, 1E5, 1E5]
terminationConditions = ['KEmin', 'maxt', 'maxt', 'maxt', 'maxt', 'maxt', 'maxt', 'maxt']
groupAnchorings = ['all', 'all', 'all', 'all', 'all', 'all', 'all', 'none']
DCDFileWrite = True
DCDFileDumpsteps = [0]

# ---=== Chromophore Parameters ===---

AARigidBodySpecies = {\
}
CGSiteSpecies = {\
CGELECTRONICSPECIES
}
maximumHopDistance = 10.0
removeORCAInputs = True
removeORCAOutputs = True
chromophoreLength = 3

# ---=== Chromophore Energy Scaling Parameters ===---

literatureHOMO = None
literatureLUMO = None
targetDoSSTDHOMO = None
targetDoSSTDLUMO = None


# ---=== Kinetic Monte Carlo Parameters ===---

systemTemperature = 290
numberOfHolesPerSimulationTime = 100000
numberOfElectronsPerSimulationTime = 0
hopLimit = 0
simulationTimes = [1.00e-12, 3.98e-12, 1.58e-11, 6.31e-11, 2.51e-10, 1.00e-9]
recordCarrierHistory = True
reorganisationEnergyDonor = None
reorganisationEnergyAcceptor = None
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
