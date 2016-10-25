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
executeZINDOS = True
executeCalculateTransferIntegrals = True
executeCalculateMobility = True

# ---=== Fine Graining Parameters ===---

CGToTemplateDirs = {\
CGTOTEMPLATEDIRS
}
CGToTemplateFiles = {\
CGTOTEMPLATEFILES
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
moleculeTerminatingUnits = {\
'H1':[0,0,0]
}
moleculeTerminatingBonds = [\
]
moleculeTerminatingConnections = [\
TERMINATINGCONNECTIONS
]

# ---=== Forcefield Parameters ===---
pairRCut = 10
pairDPDGammaVal = 0.0
# --== Lennard-Jones Pair ==--
ljCoeffs = [\
LJPAIRCOEFFS
]
# --== Disipative Particle Dynamics Pair ==--
dpdCoeffs = [\
DPDPAIRCOEFFS
]
# --== Bond ==--
bondCoeffs = [\
BONDCOEFFS
]
# --== Angle ==--
angleCoeffs = [\
ANGLECOEFFS
]
# --== Dihedral ==--
dihedralCoeffs = [\
DIHEDRALCOEFFS
]
# --== Improper ==--
improperCoeffs = [\
IMPROPERCOEFFS
]

# ---=== Molecular Dynamics Phase Parameters ===---
numberOfPhases = 8
temperatures = [1.0]
pairType = ['none', 'dpd', 'lj', 'lj', 'lj', 'lj', 'lj', 'lj']
bondType = ['harmonic']
angleType = ['harmonic']
dihedralType = ['table']
integrationTargets = ['all']
timesteps = [1E-3, 1E-3, 1E-9, 5E-9, 1E-8, 1E-7, 1E-6, 1E-5]
phaseDurations = [1E5, 1E4, 1E3, 1E3, 1E3, 1E4, 1E5, 1E5]
terminationConditions = ['KEmin', 'maxt', 'maxt', 'maxt', 'maxt', 'maxt', 'maxt', 'maxt']
groupAnchoring = ['all', 'all', 'all', 'all', 'all', 'all', 'all', 'none']


# ---=== Begin run ===---
parameterFile = __file__

if __name__ == "__main__":
    import runMorphCT

    parameterNames = [i for i in dir() if (not i.startswith('__')) and (i not in ['runMorphCT', 'os', 'sys'])]
    parameters = {}
    for name in parameterNames:
        parameters[name] = locals()[name]
    runMorphCT.simulation(**parameters) # Execute MorphCT using these simulation parameters
