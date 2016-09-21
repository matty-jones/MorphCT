import runMorphCT

# ---=== Directory and File Structure ===---
inputDir = INPUTDIR
outputDir = OUTPUTDIR


# ---=== Input Morphology Details ===---
morphology = INPUTMORPHOLOGY
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

repeatUnitTemplateDirectory = TEMPLATEDIR
repeatUnitTemplateFile = AATEMPLATEFILE
CGToTemplateAAIDs = {\
CGTOTEMPLATEAAIDS
}
CGToTemplateBonds = {\
CGTOTEMPLATEBONDS
}
### NEED TO INCLUDE RIGID BODIES HERE ###

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
numberOfPhases = 6
temperatures = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
pairType = ['none', 'dpd', 'lj', 'lj', 'lj', 'lj']
bondType = 'harmonic'
angleType = 'harmonic'
dihedralType = 'table'
integrationTargets = ['all', 'sidechains', 'all', 'all', 'all', 'all']
timesteps = [1E-3, 1E-3, 1E-9, 1E-7, 1E-6, 1E-5]
phaseDurations = [1E3, 1E4, 1E2, 1E2, 1E5, 1E5]
terminationConditions = ['KEmin', 'maxt', 'maxt', 'maxt', 'maxt', 'maxt']
groupAnchoring = ['all', 'all', 'all', 'all', 'all', 'none']


# ---=== Begin run ===---
parameterFile = __file__

if __name__ == "__main__":
    parameterNames = [i for i in dir() if (not i.startswith('__')) and (i not in ['runMorphCT'])]
    parameters = {}
    for name in parameterNames:
        parameters[name] = locals()[name]
    runMorphCT.simulation(**parameters) # Execute MorphCT using these simulation parameters
