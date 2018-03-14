# ---==============================================---
# ---======== Directory and File Structure ========---
# ---==============================================---

'''
The following parameters describe the location of the required input and
eventual output files. Note that Morph directories must be specified, even
for device files to ensure the correct device components are loaded.
'''
inputMorphDir = '/Users/mattyjones/GoogleDrive/Boise/Code/MorphCT/inputMolMorphs'
outputMorphDir = '/Users/mattyjones/GoogleDrive/Boise/Code/MorphCT/outputMolFiles'
inputDeviceDir = '/Users/mattyjones/GoogleDrive/Boise/Code/MorphCT/inputDeviceMorphs'
outputDeviceDir = '/Users/mattyjones/GoogleDrive/Boise/Code/MorphCT/outputDeviceFiles'

# ---==============================================---
# ---========== Input Morphology Details ==========---
# ---==============================================---

morphology = None
'''
The name of the morphology to run (needs to match the name of the xml (including file extension) in inputMorphDir.
Can be set to None if only running device simulations (and molecular sims are already completed)
'''

inputSigma = 1.0
'''
The sigma value to use (in Angstroems) - necessary for the Quantum Chemical Calculations.
Outputs will be unwrapped to the Angstroem length scale.
'''

deviceMorphology = 'testCrystalBilayerLog'
'''
The name of the device morphology to use (must match a directory in the inputDeviceDir).
'''

deviceComponents = {
    0: 'donorCrystal',
    1: 'acceptorCrystal',
    2: 'mixedCrystalBilayer',
}
'''
The names of the corresponding device components/moieties (must match morphologies with calculateTransferIntegrals
already executed in the outputMorphDir)
'''

overwriteCurrentData = True  # A flag that overwrites data that already exists in the relevant output directory.

# ---==============================================---
# ---============= Execution Modules ==============---
# ---==============================================---
'''
The following section allows the user to select which MorphCT modules they would like to run.
Comments after each module describe the prerequisites that must be run first.
'''
executeFinegraining = False                 # Requires: None
executeMolecularDynamics = False            # Requires: FineGraining
executeExtractMolecules = False             # Requires: MolecularDynamics
executeObtainChromophores = False           # Requires: Atomistic Morphology, or MolecularDynamics
executeZINDO = False                        # Requires: ObtainChromophores
executeCalculateTransferIntegrals = False   # Requires: ExecuteZINDO
executeCalculateMobility = False            # Requires: CalculateTransferIntegrals
executeDeviceSimulation = True              # Requires: CalculateTransferIntegrals for all deviceComponents

# ---==============================================---
# ---========== Fine Graining Parameters ==========---
# ---==============================================---

# The following parameters are used to map an input coarse-grained morphology to particular template species.
CGToTemplateDirs = {
}
'''
The directory that points to the relevant template files
(KEYS = CG site type, VALUES = directory to use)
'''

CGToTemplateFiles = {
}
'''
# The xml files to use as template
        # (KEYS = CG site type, VALUES = xml template to use)
'''

CGToTemplateForceFields = {
}
'''
The forcefield files to use as templates
(KEYS = CG site type, VALUES = forcefield xml)
'''

CGToTemplateAAIDs = {
}
'''
The mapping of coarse-grained sites to the AAIDs given in the template file.
For example, in P3HT the CG site 'A' maps to the thiophene ring which corresponds
to AAIDs 0, 1, 2, 3, 4, 24 in the template file.
(KEYS = CG site type, VALUES = list of corresponding AAIDs)
'''

CGToTemplateBonds = {
}
'''
The mapping of coarse-grained bonds to atomstic bonds. For example, in P3HT 'bondB' in the CG input morphology
corresponds to the 'C2-C3' bond between atoms #2 and #5, so the dictionary item look like this:
'bondB': ['C2-C3', 2, 5]
(KEYS = CG bond name, VALUES = ['AA bond name', AAID1, AAID2])
'''

rigidBodySites = {
}
'''
The description of the rigid bodies (if present) in the system. For example, in P3HT the thiophene
ring is rigid and so we define a rigidBodySite as:
'A': [0, 1, 2, 3, 4]. Note that any additional atoms belonging to the coarse-grained site
that are not included here will be treated as flexible.
(KEYS = CG site type, VALUES = list of AAIDs belonging to the rigid body)
'''

additionalConstraints = [
]
'''
A list of the constraints that need to be added to the system that describe the constraints between
coarse-grained sites. For example, in P3HT an additional bonded constraint is required between monomers,
between atom #3 on one monomer and atom #0 on the next monomer. Since there are 25 atoms defined in
the template file for each monomer, the entry looks like this:
['A', 'C1-C10', 3, 25]. Note that this treatment works for bonds, angles, and dihedrals and the
constraint is configured according to the length of the additionalConstraints element.
(['CG Site Name', 'Constraint name', AAID1, AAID2,...])
'''

moleculeTerminatingConnections = {
}
'''
This dictionary shows how MorphCT should terminate the chromophores if they exist at the end of the polymer chain.
This is important, because templates usually describe a `middle' monomer that exists away from the ends, which
require special treatment. The notation is identical to the addHydrogensToUA analysis script,
where {'CA': [[2, 1]]} means "for 'CA' atoms with 2 bonds, add a single hydrogen".
'''

# ---==============================================---
# ---=========== Forcefield Parameters ============---
# ---==============================================---

pairRCut = 10.0  # Cut-off for pair potentials during executeMolecularDynamics

pairDPDGammaVal = 0.0
'''
The value of the dissipative gamma to use during the DPD phase of executeMolecularDynamics (unused if no DPD)
'''

# ---==============================================---
# ---===== Molecular Dynamics Phase Parameters ====---
# ---==============================================---
'''
The following parameters describe the structure of the executeMolecularDynamics phase. As a general rule, the
length of each parameter should be equal to the specified numberOfPhases, however if only one element is
specified then that value will be used for all phases.
'''

numberOfPhases = 8              # The number of MD phases to run
temperatures = [1.0]            # The system temperature (dimensionless
taus = [1.0]                    # The thermostat coupling
pairTypes = ['none', 'dpd', 'lj', 'lj', 'lj', 'lj', 'lj', 'lj']
'''
The pair interactions to use (permitted: 'none', 'dpd', 'lj')
'''

bondTypes = ['harmonic']        # The bond constraint equation to use
angleTypes = ['harmonic']       # The angle constraint equation to use
dihedralTypes = ['opls']        # The dihedral constraint equation to use
integrationTargets = ['all']    # The atoms to be integrated (permitted: 'all', <specific CG site name>)
timesteps = [1E-3, 1E-3, 1E-10, 1E-9, 1E-8, 1E-7, 1E-6, 1E-5]
'''
Timestep values for each phase (in dimensionless time units)
'''

durations = [1E5, 1E4, 1E3, 1E3, 1E3, 1E4, 1E5, 1E5]  # Simulation durations for each phase (in # of timesteps)
terminationConditions = ['KEmin', 'maxt', 'maxt', 'maxt', 'maxt', 'maxt', 'maxt', 'maxt']
'''
The conditions under which the simulation finishes. Incorporated to allow simulations to finish when
they hit the minimum KE. (permitted: 'KEmin, maxt')
'''

groupAnchorings = ['all', 'all', 'all', 'all', 'all', 'all', 'all', 'none']
'''
Usually, when 'all' is selected, the rigid body CoMs are strongly bonded with 0 equilibration distance to
the original coarse-grained site positions. This constraint can be relaxed by specifying either particular
CG sites to constrain (e.g. 'A,B') or 'none'
'''

DCDFileWrite = True             # Select whether to write a DCD file
DCDFileDumpsteps = [0]
'''
Select the frequency of DCD writes. [0] defaults to current phase duration / 100 (i.e. 100 total trajectory frames).
Otherwise the value corresponds to the DCD dump period.
'''

# ---==============================================---
# ---============ Chromophore Parameters ==========---
# ---==============================================---
'''
For AARigidBodySpecies and CGSiteSpecies there are 3 modes of operation:

   1) len(AARigidBodySpecies) = 0, len(CGSiteSpecies) > 1:
      This is normal operation, where the CGSiteSpecies is a dictionary
      that maps the <CGSiteType>: <ElectronicType>. This requires the
      fine-graining module to have been run by MorphCT. For example in
      coarse-grained P3HT, CGSiteSpecies = {'A': 'Donor', 'B': 'None', 'C': 'None'}.

   2) len(CGSiteSpecies) == 1:
      This is the operation for a morphology where there is only a single
      type of electronic species (e.g. neat small molecule system). All
      chromophores will be set to this species type and the key does not matter.
      For example in neat Perylene, CGSiteSpecies = {'A': Donor}

   3) len(AARigidBodySpecies) > 0, len(CGSiteSpecies) == 1:
      This is the operation for a more complex atomistic morphology, where multiple species
      may be present (e.g. block copolymer, or small molecule blend). MorphCT uses the HOOMD
      rigid bodies to decide which chromophores are which species. Note that non-numpy
      ranges work too. For example,
      AARigidBodySpecies = {'Donor': range(0,100), 'Acceptor': range(100,200)}
'''

AARigidBodySpecies = {
}
CGSiteSpecies = {
}
useVoronoiNeighbours = True
'''
If true, the voronoi analysis will be used to find hopping neighbours.
Otherwise, defaults to cut-off with distances specified below.
'''

maximumHoleHopDistance = 10.0       # Set maximum hop distance for holes (unused if useVoronoiNeighbours == True)
maximumElectronHopDistance = 10.0   # Set maximum hop distance for electrons (unused if useVoronoiNeighbours == True)
permitHopsThroughOpposingChromophores = False
'''
If true, obtainChromophores will ignore the opposing chromophore type when determining the neighbours of a chromophore.
This manifests as the ability of a carrier to hop from one donor to another 'through' an acceptor (or vice versa) as
if the intermediate chromophore was not there.
'''

removeORCAInputs = True
'''
Set ORCA input files to be deleted after the transfer integrals have been correctly obtained (recommended)
'''

removeORCAOutputs = True
'''
Set ORCA output files to be deleted after the transfer integrals have been correctly obtained (not
recommended for systems with < 10,000 chromophore outputs as it makes it easier to debug)
'''

chromophoreLength = 3
'''
This parameter is no longer supported on the master branch due to lack of testing.
Check the MorphCT "variableChromoLengths" if you'd like to.
'''

# ---=== Chromophore Energy Scaling Parameters ===---
'''
targetDoSSTD (float, eV)
    Gaussian width of the donor HOMO/LUMO density of states (100 meV recommended for polymers,
    small molecules can be less)

VRHDelocalisation (float, meters)
    The carrier delocalisation length for VRH (in m) such that \alpha -> \alpha * (1.0 / VRHDelocalisation)

literatureMO (float, eV)
    Literature value of the HOMO/LUMO of the donor/acceptor

species (string)
    Donor type or Acceptor type

reorganisationEnergy (float, eV)
'''


chromophore_species = {
  "D1": {
    "literatureMO": -5.0,
    "targetDoSSTD": 0.1,
    "reorganisationEnergy": 0.3064,
    "species": "Donor",
    "VRHDelocalisation": 2e-10,
  },
  "D2": {
    "literatureMO": -5.0,
    "targetDoSSTD": 0.1,
    "reorganisationEnergy": 0.3064,
    "species": "Donor",
    "VRHDelocalisation": 2e-10,
  },
  "A1": {
    "literatureMO": -4.7,
    "targetDoSSTD": 0.1,
    "reorganisationEnergy": 0.1496,
    "species": "Acceptor",
    "VRHDelocalisation": 4e-10,
  },
  "A2": {
    "literatureMO": -3.7,
    "targetDoSSTD": 0.1,
    "reorganisationEnergy": 0.2496,
    "species": "Acceptor",
    "VRHDelocalisation": 4e-10,
  },
}

useKoopmansApproximation = False
'''
Treat all chromophores as identical (i.e. same HOMO and LUMO levels). DeltaEij will be set to zero both in the transfer
integral calculation and in the Marcus hopping rate. Mobilities will be significantly higher than experiment, but it
will permit the user to only consider the effect of morphology on charge transport.
'''

koopmansHoppingPrefactor = 1E-3
'''
When Koopmans' approximation is active, this hopping prefactor (which is dependent on energetic disorder) will be
applied, which can be tuned to obtain better agreement with experimental mobilities. NOTE: This hopping prefactor
is only considered in the mobility simulations. For device simulations, only the hoppingPrefactor parameter is
used (look in the Device KMC Parameters section).
'''

# ---==============================================---
# ---=== General Kinetic Monte Carlo Parameters ===---
# ---==============================================---

# ---=== Universal KMC Parameters ===---
# The following parameters are universally relevant for both morphology and device simulations
systemTemperature = 290                 # Device temperature (RTP)
useSimpleEnergeticPenalty = False
'''
Replaces the exponential term in the Marcus hopping rate with a simple Boltzmann penalty for hops upstream in
energy (similar to that in Miller Abraham's hopping). Mathematically:
\exp^{\frac{- ((\Delta E_{ij} + \lambda_{ij})**2}{4 \lambda_{ij} k_{B} T}}
becomes
\exp^{\frac{- \Delta E_{ij}}{k_{B} T}}.
'''

recordCarrierHistory = True
'''
Required to plot connectivity graphs, but drastically increases pickle size, processing time and memory usage
(although uses scipy's sparse matrices so it's not too bad). The following parameters are concerned with the
presence of 'variable range hopping', i.e. providing an energetic penalty for long hops. Variable Range Hopping
looks like this:
kij_{VRH} = vij_{NoVRH} * \exp^{- \alpha r_{ij}}
'''

useVRH = True
'''
Include an explicit consideration of hop distance into the hopping rate equation (otherwise long-range hop
decay is entirely determined by the transfer integral)
'''


# ---=== Mobility Specific KMC Parameters ===---
numberOfHolesPerSimulationTime = 0      # The total number of holes to simulate
numberOfElectronsPerSimulationTime = 0  # The total number of electrons to simulations
hopLimit = 0
'''
For testing, rather than running for a specific simulation time, terminate the carrier after this number of hops
'''

simulationTimes = [1.00e-12, 3.98e-12, 1.58e-11, 6.31e-11, 2.51e-10, 1.00e-9]
'''
The termination condition for the mobility KMC simulations
'''

combineKMCResults = True
'''
This flag combines the KMC results from multiple cores before termination. Often this doesn't finish in time
due to cluster configuration, but there's an anlysis script that does it for you (combineKMC) if it doesn't work.
'''

useAverageHopRates = False
averageIntraHopRate = 8.07E14
averageInterHopRate = 3.92E14

# ---=== Device Kinetic Monte Carlo Parameters ===---
# The following parameters are specific to the device simulations only
#           Material
absorptionCoefficient = 1.3E4   # cm^{-1}
relativePermittivity = 3
donorHOMO = -5.3                # eV
acceptorLUMO = -3.9             # eV
cathodeWorkFunction = -4.2      # eV
anodeWorkFunction = -5.0        # eV
recombinationRate = 1E9         # s^{-1}
coulombCaptureRadius = 1E-9     # m (carriers will recombine if they are within this distance)
wrapDeviceXY = False
'''
Carriers leaving the device through the X/Y axes will wrap back onto the other side. Otherwise these hops are forbidden.
'''

# Charged Species
excitonLifetime = 0.5E-9        # s
forsterRadius = 4.3E-9          # m
hoppingPrefactor = 1E-4         # Prefactor to slow down carrier and exciton hops
MAPrefactor = 1E11  # Prefactor to calibrate the Miller Abrahams dark-current injection hops (1E-6 * mean of rate dist)
MALocalisationRadius = 1E-9     # m
# External Parameters
electricalField = 1E7           # V/m
incidentFlux = 10000            # mW/cm^{2}
incidentWavelength = 500E-9     # m
# Simulation
voltageSweep = [-0.1, 0.4, 0.8]  # Voltage values to sweep over for the generation of the J-V curve.
morphologyCellSize = 1.00E-8    # m (This value comes from the average of the moietyTypes used for testing)
minimumNumberOfPhotoinjections = 100  # Determines the end of the device simulations
fastestEventAllowed = 1E-15     # s, Events with tau < x will be discarded, can be None to use defaults
slowestEventAllowed = 1E-6      # s, Events with tau > x will be discarded, can be None to use defaults
disableDarkInjection = True     # Disables carrier injection from the anode/cathode
disableCoulombic = True         # Stops carriers `feeling' each other through Coulombic interactions
#           Logging
outputLogToSTDOUT = True
'''
Divert log output to the terminal rather than the KMClog files (useful for testing or when running only one voltage
on a single core)
'''

# ---==============================================---
# ---================= Begin run ==================---
# ---==============================================---

parameterFile = __file__

if __name__ == "__main__":
    from morphct import runMorphCT
    from morphct.code import helperFunctions

    procIDs = helperFunctions.getCPUCores()
    parameterNames = [i for i in dir() if (not i.startswith('__')) and
                      (i not in ['runMorphCT', 'helperFunctions', 'sys'])]
    parameters = {}
    for name in parameterNames:
        parameters[name] = locals()[name]
    runMorphCT.simulation(**parameters)  # Execute MorphCT using these simulation parameters
