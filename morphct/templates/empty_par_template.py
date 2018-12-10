import os as __os

# ---==============================================---
# ---======== Directory and File Structure ========---
# ---==============================================---

"""
The following parameters describe the location of the required input and
eventual output files. Note that morph directories must be specified, even
for device files to ensure the correct device components are loaded.
A seperate orca directory can be specified for the QCCs. It is recommended
that /dev/shm or similar shared-memory directory be used to dramatically
improve the speed of orca calculations. Otherwise, setting the output_orca_dir
to None will use the output_morph_dir instead.
"""
input_morph_dir =
output_morph_dir =
output_orca_dir =
input_device_dir =
output_device_dir =

# ---==============================================---
# ---========== Input Morphology Details ==========---
# ---==============================================---
"""
The name of the morphology to run (needs to match the name of the xml (including file extension) in input_morph_dir.
Can be set to None if only running device simulations (and molecular sims are already completed)
"""
morphology =

"""
The sigma value to use (in Angstroems) - necessary for the quantum chemical calculations.
Outputs will be unwrapped to the angstroem length scale.
"""
input_sigma =

"""
The name of the device morphology to use (must match a directory in the input_device_dir).
"""
device_morphology =

"""
The names of the corresponding device components/moieties (must match morphologies with calculate_transfer_integrals
already executed in the output_morph_dir)
"""
device_components =
"""
The boolean that ignores preveiously calculated data for this morphology, and instead overwrites it in the relevant
output directory.
"""
overwrite_current_data =

"""
An integer to set the master random seed in definitions.py. From this seed, all child process seeds will be spawned.
Default is None, which will generate a seed based on the system time at runtime. Setting the seed to a 32-bit integer
will ensure that runs with identical parameter files will generate identical data.
"""
random_seed_override =


# ---==============================================---
# ---============= Execution Modules ==============---
# ---==============================================---
"""
The following section allows the user to select which MorphCT modules they would like to run.
comments after each module describe the prerequisites that must be run first.
"""
execute_fine_graining =
execute_molecular_dynamics =
execute_obtain_chromophores =
execute_ZINDO =
execute_calculate_transfer_integrals =
execute_calculate_mobility =
execute_device_simulation =

# ---==============================================---
# ---========== Fine Graining Parameters ==========---
# ---==============================================---
"""
The following section allows the user to map an input coarse-grained morphology to a particular template species
"""

"""
The directory that points to the relevant template files
(KEYS = CG site type, VALUES = directory to use)
"""
CG_to_template_dirs =

"""
The xml files to use as template
(KEYS = CG site type, VALUES = xml template to use)
"""
CG_to_template_files =

"""
The forcefield files to use as templates
(KEYS = CG site type, VALUES = forcefield xml)
"""
CG_to_template_force_fields =

"""
The mapping of coarse-grained sites to the AAIDs given in the template file.
For example, in P3HT the CG site 'A' maps to the thiophene ring which corresponds
to AAIDs 0, 1, 2, 3, 4, 24 in the template file.
(KEYS = CG site type, VALUES = list of corresponding AAIDs)
"""
CG_to_template_AAIDs =

"""
The mapping of coarse-grained bonds to atomstic bonds. For example, in P3HT 'bondB' in the CG input morphology
corresponds to the 'C2-C3' bond between atoms #2 and #5, so the dictionary item look like this:
'bondB': ['C2-C3', 2, 5]
(KEYS = CG bond name, VALUES = ['AA bond name', AAID1, AAID2])
"""
CG_to_template_bonds =

"""
The description of the rigid bodies (if present) in the system. For example, in P3HT the thiophene
ring is rigid and so we define a rigid_body_site as:
'A': [0, 1, 2, 3, 4]. Note that any additional atoms belonging to the coarse-grained site
that are not included here will be treated as flexible.
(KEYS = CG site type, VALUES = list of AAIDs belonging to the rigid body)
"""
rigid_body_sites =

"""
A list of the constraints that need to be added to the system that describe the constraints between
coarse-grained sites. For example, in P3HT an additional bonded constraint is required between monomers,
between atom #3 on one monomer and atom #0 on the next monomer. Since there are 25 atoms defined in
the template file for each monomer, the entry looks like this:
['A', 'C1-C10', 3, 25]. Note that this treatment works for bonds, angles, and dihedrals and the
constraint is configured according to the length of the additional_constraints element.
(['CG Site Name', 'Constraint Name', AAID1, AAID2,...])
"""
additional_constraints =

"""
This dictionary shows how MorphCT should terminate the chromophores if they exist at the end of the polymer chain.
This is important, because templates usually describe a `middle' monomer that exists away from the ends, which
require special treatment. The notation is identical to the add_hydrogens_to_ua analysis script,
where {'CA': [[2, 1]]} means "for 'CA' atoms with 2 bonds, add a single hydrogen".
"""
molecule_terminating_connections =

# ---==============================================---
# ---=========== Forcefield Parameters ============---
# ---==============================================---
"""
The following section allows the user to calibrate some parameters for the pair potentials
"""

"""
The cut-off value for pair potentials during execute_molecular_dynamics
"""
pair_r_cut =

"""
The value of the dissipative gamma to use during the DPD phase of execute_molecular_dynamics (unused if no DPD)
"""
pair_dpd_gamma_val =

# ---==============================================---
# ---===== Molecular Dynamics Phase Parameters ====---
# ---==============================================---
"""
The following parameters describe the structure of the execute_molecular_dynamics phase. As a general rule, the
length of each parameter should be equal to the specified number_of_phases, however if only one element is
specified then that value will be used for all phases.
"""

"""
The number of MD phases to run
"""
number_of_phases =

"""
The dimensionless system temperature to perform MD at
"""
temperatures =

"""
The dimensionless thermostat coupling
"""
taus =

"""
The pair interactions to use (permitted: 'none', 'dpd', 'lj')
"""
pair_types =

"""
The bond constraints to use (permitted: 'harmonic')
"""
bond_types =

"""
The angle constraints to use (permitted: 'harmonic')
"""
angle_types =

"""
The dihedral constraints to use (permitted: 'opls', 'table')
"""
dihedral_types =

"""
The atoms to be included in the integration steps (permitted: 'all', <specific CG site name>)
"""
integration_targets =

"""
Timestep values for each phase (in dimensionless time units)
"""
timesteps =

"""
Simulation durations for each phase (in # of timesteps)
"""
durations =

"""
The conditions under which the simulation finishes. Incorporated to allow simulations to finish when
they hit the minimum KE. (permitted: 'ke_min, max_t')
"""
termination_conditions =

"""
Usually, when 'all' is selected, the rigid body COMs are strongly bonded with 0 equilibration distance to
the original coarse-grained site positions. This constraint can be relaxed by specifying either particular
CG sites to constrain (e.g. 'A,B') or 'none'
"""
group_anchorings =

"""
Boolean to decide whether to write a DCD file
"""
dcd_file_write =

"""
Select the frequency of DCD writes. [0] defaults to current phase duration / 100 (i.e. 100 total trajectory frames).
Otherwise the value corresponds to the DCD dump period.
"""
dcd_file_dumpsteps =

# ---==============================================---
# ---============ Chromophore Parameters ==========---
# ---==============================================---
"""
For AA_rigid_body_species and CG_site_species there are 3 modes of operation:

   1) len(AA_rigid_body_species) = 0, len(CG_site_species) > 1:
      This is normal operation, where the CG_site_species is a dictionary
      that maps the <CGSiteType>: <electronicType>. this requires the
      fine-graining module to have been run by MorphCT. For example in
      coarse-grained P3HT, CG_site_species = {'A': 'donor', 'B': 'none', 'C': 'none'}.

   2) len(AA_rigid_body_species) = 0, len(CG_site_species) == 1:
      This is the operation for a morphology where there is only a single
      type of electronic species (e.g. neat small molecule system). All
      chromophores will be set to this species type and the key does not matter.
      For example in neat perylene, CG_site_species = {'A': donor}

   3) len(AA_rigid_body_species) > 0:
      This is the operation for a more complex atomistic morphology, where multiple species
      may be present (e.g. block copolymer, or small molecule blend). MorphCT uses the HOOMD
      rigid bodies to decide which chromophores are which species. Note that non-numpy
      ranges work too. For example,
      AA_rigid_body_species = {'donor': range(0,100), 'acceptor': range(100,200)}
"""
AA_rigid_body_species =
CG_site_species =

"""
If true, the voronoi analysis will be used to find hopping neighbours.
Otherwise, defaults to cut-off with distances specified below.
"""
use_voronoi_neighbours =

"""
If use_voronoi_neighbours is False, the maximum hop distance for holes
"""
maximum_hole_hop_distance =

"""
If use_voronoi_neighbours is False, the maximum hop distance for electrons
"""
maximum_electron_hop_distance =

"""
If True, obtain_chromophores will ignore the opposing chromophore type when determining the neighbours of a
chromophore. This manifests as the ability of a carrier to hop from one donor to another 'through' an acceptor (or
vice versa) as if the intermediate chromophore was not there.
"""
permit_hops_through_opposing_chromophores =

"""
Set ORCA input files to be deleted after the transfer integrals have been correctly obtained (recommended)
"""
remove_orca_inputs =

"""
Set ORCA output files to be deleted after the transfer integrals have been correctly obtained (not
recommended for systems with < 10,000 chromophore outputs as it makes it easier to debug)
"""
remove_orca_outputs =

"""
This parameter is no longer supported on the master branch due to lack of testing.
Check the MorphCT "variable_chromo_lengths" if you'd like to.
"""
chromophore_length =

# ---==============================================---
# ---=== Chromophore Energy Scaling Parameters ====---
# ---==============================================---
"""
target_DOS_std (float, e_V)
    Gaussian width of the donor HOMO/LUMO density of states (100 meV recommended for polymers,
    small molecules can be less)

VRH_delocalisation (float, meters)
    The carrier delocalisation length for VRH (in m) such that \alpha -> \alpha * (1.0 / VRH_delocalisation)

literature_MO (float, e_V)
    Literature value of the HOMO/LUMO of the donor/acceptor

species (string)
    donor type or acceptor type

reorganisation_energy (float, e_V)
"""
chromophore_species =

"""
Treat all chromophores as identical (i.e. same HOMO and LUMO levels). Delta_Eij will be set to zero both in the transfer
integral calculation and in the marcus hopping rate. mobilities will be significantly higher than experiment, but it
will permit the user to only consider the effect of morphology on charge transport.
"""
use_koopmans_approximation =

"""
When Koopmans' approximation is active, this hopping prefactor (which is dependent on energetic disorder) will be
applied, which can be tuned to obtain better agreement with experimental mobilities. Note: This hopping prefactor
is only considered in the mobility simulations. For device simulations, only the hopping_prefactor parameter is
used (look in the device KMC parameters section).
"""
koopmans_hopping_prefactor =

# ---==============================================---
# ---=== General Kinetic Monte Carlo Parameters ===---
# ---==============================================---

# ---=== Universal KMC Parameters ===---
"""
The following parameters are universally relevant for both mobility and device simulations
"""

"""
The device temperature when the KMC simulations are performed
"""
system_temperature =

"""
Replaces the exponential term in the Marcus hopping rate with a simple Boltzmann penalty for hops upstream in
energy (similar to that in Miller Abraham's hopping). Mathematically:
\exp^{\frac{- ((\Delta E_{ij} + \lambda_{ij})**2}{4 \lambda_{ij} k_{B} T}}
becomes
\exp^{\frac{- \Delta E_{ij}}{k_{B} T}}.
"""
use_simple_energetic_penalty =

"""
Required to plot connectivity graphs, but drastically increases pickle size, processing time and memory usage
(although uses scipy's sparse matrices so it's not too bad). The following parameters are concerned with the
presence of 'variable range hopping', i.e. providing an energetic penalty for long hops. Variable range hopping
looks like this:
kij_{VRH} = vij_{NoVRH} * \exp^{- \alpha r_{ij}}
"""
record_carrier_history =

"""
Include an explicit consideration of hop distance into the hopping rate equation (otherwise long-range hop
decay is entirely determined by the transfer integral)
"""
use_VRH =


# ---=== Mobility Specific KMC Parameters ===---
"""
The following parameters are relevant only for the mobility KMC simulations
"""

"""
The total number of holes to simulate
"""
number_of_holes_per_simulation_time =

"""
The total number of holes to simulate
"""
number_of_electrons_per_simulation_time =

"""
For testing, rather than running for a specific simulation time, terminate the carrier after this number of hops
"""
hop_limit =

"""
The termination condition for the mobility KMC simulations
"""
simulation_times =

"""
This flag combines the KMC results from multiple cores before termination. Often this doesn't finish in time
due to cluster configuration, but there's an anlysis script that does it for you (combine_KMC) if it doesn't work.
"""
combine_KMC_results =

"""
This flag allows the user to define a fixed intra- and inter-molecular hopping rate (sometimes useful when explore
mobility landscape). If true, all carriers will have the hop rates defined below
"""
use_average_hop_rates =

"""
The fixed intra-molecular hop rate to use when use_average_hop_rates is True
"""
average_intra_hop_rate =

"""
The fixed inter-molecular hop rate to use when use_average_hop_rates is True
"""
average_inter_hop_rate =

# ---=== Device Kinetic Monte Carlo Parameters ===---
"""
The following parameters are relevant to the device KMC simulations only
"""

#           Material
"""
The following parameters describe the simulated material
"""

"""
The absorption coefficient in units of cm^{-1}
"""
absorption_coefficient =

"""
The relative permittivity of the material
"""
relative_permittivity =

"""
The HOMO level of the simulated donor material in units of eV
"""
donor_HOMO =

"""
The LUMO level of the simulated acceptor materials in units of eV
"""
acceptor_LUMO =

"""
The workfunction of the cathode in units of eV
"""
cathode_work_function =

"""
The workfunction of the anode in units of eV
"""
anode_work_function =

"""
The recombination rate of opposing carriers in units of s^{-1}
"""
recombination_rate =

"""
The coulomb capture radius of carriers in units of m. A recombination time will be calculated between two opposing
carriers if they are separated by less than this distance. If the carriers are still within this radius when the
recombination event takes place, the carriers recombine and are removed from the simulations.
"""
coulomb_capture_radius =

"""
Carriers leaving the device through the X/Y axes will wrap back onto the other side. otherwise these hops are forbidden.
"""
wrap_device_xy =

# Charged Species
"""
The following parameters describe the charged species
"""

"""
The exciton lifetime in units of s
"""
exciton_lifetime =

"""
The Forster radius for the exciton in units of m
"""
forster_radius =

"""
A multiplicative hopping prefactor that slows down all carrier and exciton hops to bring their rates more in keeping
with the other event types (otherwise the order of magnitude time discrepancy is way too high and nothing actually
happens in the device).
"""
hopping_prefactor =

"""
The Miller Abrahams prefactor that controls the rate of dark-current injection hops.
"""
MA_prefactor =

"""
The Miller Abrahams localisation radius in units of m.
"""
MA_localisation_radius =

# External Parameters
"""
The following parameters describe the incident photons and other factors external to the device
"""

"""
The incident flux in units of m_W/cm^{2}
"""
incident_flux =

"""
The wavelength of the incident light in units of m
"""
incident_wavelength =

# Simulation
"""
The following parameters describe the remainder of the simulation
"""

"""
The voltage values to sweep over for the generation of the J-V curve
"""
voltage_sweep =

"""
The average size of the constituent molecular morphologies within the device (used to calculate separations over
multiple moieties) in units of m
"""
morphology_cell_size =

"""
The total number of photoinjections to simulate. After this, the device simulation will end
"""
minimum_number_of_photoinjections =

"""
The minimum event timescale to permit in the system. Events with tau < fastest_event_allowed will be discarded. Can be
set to None to use the defaults
"""
fastest_event_allowed =

"""
The maximum event timescale to permit in the system. Events with tau > slowest_event_allowed will be discarded. Can be
set to None to use the defaults
"""
slowest_event_allowed =

"""
Boolean that can be used to toggle carrier injection from the anode/cathode
"""
disable_dark_injection =

"""
Boolean that can be used to toggle Coulombic interactions between carriers in the device at the same time
"""
disable_coulombic =

"""
Divert log output to the terminal rather than the KMC_log files (useful for testing or when running only one voltage
on a single core)
"""
output_log_to_stdout =

# ---==============================================---
# ---================= Begin run ==================---
# ---==============================================---


if __name__ == "__main__":
    from morphct import run_MorphCT
    from morphct.code import helper_functions as hf

    parameter_file = __os.path.realpath(__file__)
    proc_IDs = hf.get_CPU_cores()
    parameter_names = [
        i
        for i in dir()
        if (not i.startswith("__"))
        and (not i.startswith("@"))
        and (i not in ["run_MorphCT", "helper_functions", "hf"])
    ]
    parameters = {}
    for name in parameter_names:
        parameters[name] = locals()[name]
    run_MorphCT.simulation(
        **parameters
    )  # Execute MorphCT using these simulation parameters
