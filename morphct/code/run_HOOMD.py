import os
import sys
import numpy as np

try:
    from hoomd_script import *
except ImportError:
    print("HOOMD 1.3 NOT FOUND. FINE_GRAINING AND RUN_HOOMD WILL NOT WORK.")
from morphct.code import helper_functions as hf


class ExitHOOMD(Exception):
    """This class is raised to terminate a HOOMD simulation mid-run for
    a particular reason (e.g. minimum KE found)"""

    def __init__(self, string):
        self.string = "".join([string, " at timestep = {:d}".format(get_step())])

    def __str__(self):
        return self.string


class md_phase:
    def __init__(
        self,
        AA_morphology_dict,
        CG_morphology_dict,
        CG_to_AAID_master,
        parameter_dict,
        phase_number,
        input_file,
        output_file,
        s_scale,
        e_scale,
    ):
        """
        The MDPhase class respresents a single MD simulation using the given parameters
        """
        self.AA_morphology_dict = AA_morphology_dict
        self.CG_morphology_dict = CG_morphology_dict
        self.CG_to_AAID_master = CG_to_AAID_master
        self.input_file = input_file
        self.output_file = output_file
        self.s_scale = s_scale
        self.e_scale = e_scale
        self.phase_number = phase_number
        # Obtain the parXX.py parameters
        for key in list(parameter_dict.keys()):
            self.__dict__[key] = parameter_dict[key]
        # Get the phase-specific simulation parameters
        for key in [
            "temperatures",
            "taus",
            "pair_types",
            "bond_types",
            "angle_types",
            "dihedral_types",
            "integration_targets",
            "timesteps",
            "durations",
            "termination_conditions",
            "group_anchorings",
            "dcd_file_dumpsteps",
        ]:
            # If the phase-specific parameter is not defined for this phase
            # number then use the first one.
            if self.phase_number + 1 > len(parameter_dict[key]):
                self.__dict__[key[:-1]] = parameter_dict[key][0]
            # If the phase-specific parameter is specified for this phase then
            # use this parameter
            else:
                self.__dict__[key[:-1]] = parameter_dict[key][phase_number]
        self.system = self.load_system(input_file)
        # Determine the required groups so we can use the correct integrator
        # each time
        self.rigid_group, self.non_rigid_group, self.integration_types = (
            self.get_integration_groups()
        )
        self.output_log_file_name = os.path.join(
            self.output_morph_dir,
            os.path.splitext(self.morphology)[0],
            "morphology",
            "".join(["energies_", os.path.splitext(self.morphology)[0], ".log"]),
        )
        # Determine which quantities should be logged during the simulation
        # phase
        self.log_quantities = [
            "temperature",
            "pressure",
            "volume",
            "potential_energy",
            "kinetic_energy",
        ]
        # Set the bond coefficients
        self.get_FF_coeffs()

    def load_system(self, input_file):
        # Load the previous phases' xml for continuation
        system_xml = init.read_xml(filename=input_file)
        # A snapshot is needed in order to update the velocities
        snapshot = system_xml.take_snapshot()
        # Assign the required velocities based on the requested temperature
        updated_snapshot = self.initialize_velocities(snapshot)
        # Finally, restore the snapshot
        system_xml.restore_snapshot(updated_snapshot)
        return system_xml

    def initialize_velocities(self, snapshot):
        v = np.random.random((len(snapshot.particles.velocity), 3))
        v -= 0.5
        meanv = np.mean(v, 0)
        meanv2 = np.mean(v ** 2, 0)
        fs = np.sqrt(self.temperature / meanv2)
        # Shift the velocities such that the average is zero
        v = v - meanv
        # Scale the velocities to match the required temperature
        v *= fs
        # Assign the velocities for this MD phase
        snapshot.particles.velocity[:] = v[:]
        return snapshot

    def optimise_structure(self):
        # Activate the dumping of the trajectory dcd file
        if self.dcd_file_write is True:
            if self.dcd_file_dumpstep is not 0:
                print("Setting dcd dump step to", self.dcd_file_dumpstep)
                self.dump_dcd = dump.dcd(
                    filename=self.output_file.replace("xml", "dcd"),
                    period=self.dcd_file_dumpstep,
                    overwrite=True,
                )
            else:
                self.dump_dcd = dump.dcd(
                    filename=self.output_file.replace("xml", "dcd"),
                    period=self.duration / 100.0,
                    overwrite=True,
                )
        else:
            self.dump_dcd = None
        # Set the integrators, groups and timestep
        self.step = integrate.mode_standard(dt=self.timestep)
        self.rigid_int = integrate.nvt_rigid(
            group=self.rigid_group, T=self.temperature, tau=self.tau
        )
        self.non_rigid_int = integrate.nvt(
            group=self.non_rigid_group, T=self.temperature, tau=self.tau
        )
        # Overwrite the log file if this is the first phase, otherwise append to
        # the previous log
        if self.phase_number == 0:
            log_overwrite = True
        else:
            log_overwrite = False
        self.energy_log = analyze.log(
            filename=self.output_log_file_name,
            quantities=self.log_quantities,
            period=self.duration / 1000.0,
            overwrite=log_overwrite,
        )
        callback = None
        # Set up the callback function if the termination condition is not max_t
        if self.termination_condition.lower() == "ke_min":
            self.load_from_snapshot = False
            self.lowest_KE = 9e999
            self.KE_increased = 0
            self.first_KE_value = True
            callback = analyze.callback(
                callback=self.check_KE, period=self.duration / 1000.0
            )
        print(
            "---=== BEGINNING MOLECULAR DYNAMICS PHASE", self.phase_number + 1, "===---"
        )
        # Run the MD simulation
        try:
            run(self.duration)
        except ExitHOOMD as exitMessage:
            print(exitMessage)
        # Load the snapshot if required
        if self.termination_condition.lower() == "ke_min":
            if self.load_from_snapshot is True:
                print("Loading from snapshot...")
                self.system.restore_snapshot(self.snapshot_to_load)
                del self.snapshot_to_load
        # Create the output xml file
        self.dump_xml = dump.xml(
            filename=self.output_file,
            position=True,
            image=True,
            type=True,
            mass=True,
            diameter=True,
            body=True,
            charge=True,
            bond=True,
            angle=True,
            dihedral=True,
            improper=True,
        )
        # Clean up all references to this simulation
        self.rigid_int.disable()
        self.non_rigid_int.disable()
        del self.system, self.dump_dcd, self.step, self.rigid_int, self.non_rigid_int
        del self.rigid_group, self.non_rigid_group, callback, self.dump_xml
        del self.pair_class, self.bond_class, self.angle_class, self.dihedral_class
        del self.improper_class, self.energy_log
        init.reset()

    def check_KE(self, timestep_number):
        # Query the current kinetic energy of the system through the energy_log
        current_KE = self.energy_log.query("kinetic_energy")
        # For second and subsequent steps
        if self.first_KE_value is False:
            # Check if the current KE is greater than the minimum so far
            if current_KE >= self.lowest_KE:
                if self.KE_increased == 5:
                    # Found the lowest KE point for at least 5 dumpsteps
                    del self.first_KE_value, self.lowest_KE, self.KE_increased
                    raise ExitHOOMD("Lowest energy condition met")
                # Increment a counter that indicates how many times the KE has
                # increased since the minimum
                self.KE_increased += 1
            else:
                # At at least local KE minimum, so store snapshot
                self.KE_increased = 0
                self.load_from_snapshot = True
                self.snapshot_to_load = self.system.take_snapshot(all=True)
                self.lowest_KE = current_KE
        else:
            # Skip the first check because the KE fluctuates wildly within the
            # first dump step
            self.first_KE_value = False
        return 0

    def get_FF_coeffs(self):
        # First find all of the forcefields specified in the par file
        all_FF_names = {}
        for CG_site, directory in self.CG_to_template_dirs.items():
            FF_loc = os.path.join(directory, self.CG_to_template_force_fields[CG_site])
            if FF_loc not in list(all_FF_names.values()):
                all_FF_names[CG_site] = FF_loc
        FF_list = []
        # Then load in all of the FFs with the appropriate mappings
        for CG_site in list(all_FF_names.keys()):
            FF_list.append(
                hf.load_FF_xml(
                    all_FF_names[CG_site], mapping=self.new_type_mappings[CG_site]
                )
            )
        # Combine all of the individual, mapped FFs into one master field
        master_FF = {}
        for FF in FF_list:
            for FF_type in list(FF.keys()):
                if FF_type not in list(master_FF.keys()):
                    master_FF[FF_type] = FF[FF_type]
                else:
                    master_FF[FF_type] += FF[FF_type]
        # Finally, assign the expected variables to each value in the masterFF
        self.lj_coeffs = master_FF["lj"]
        self.dpd_coeffs = master_FF["dpd"]
        self.bond_coeffs = master_FF["bond"]
        self.angle_coeffs = master_FF["angle"]
        self.dihedral_coeffs = master_FF["dihedral"]
        self.improper_coeffs = master_FF["improper"]
        # Set Pair Coeffs
        self.pair_class = None
        if self.pair_type.lower() != "none":
            # Log the correct pairType energy
            self.log_quantities.append("".join(["pair_", self.pair_type, "_energy"]))
            # HOOMD crashes if you don't specify all pair combinations, so need
            # to make sure we do this.
            atom_types = sorted(
                list(set(self.AA_morphology_dict["type"])),
                key=lambda x: hf.convert_string_to_int(x),
            )
            all_pair_types = []
            # Create a list of all of the pairTypes to ensure that the required
            # coefficients are set
            for atom_type1 in atom_types:
                for atom_type2 in atom_types:
                    pair_type = "{0}-{1}".format(atom_type1, atom_type2)
                    reverse_pair_type = "{1}-{0}".format(atom_type1, atom_type2)
                    if (pair_type not in all_pair_types) and (
                        reverse_pair_type not in all_pair_types
                    ):
                        all_pair_types.append(pair_type)
            # Read in the pairTypes, parameters and coefficients and set them
            # for HOOMD
            if self.pair_type.lower() == "dpd":
                self.pair_class = pair.dpd(
                    r_cut=self.pair_r_cut * self.s_scale, T=self.temperature
                )
                # Use the geometric mixing rule for all possible combinations of
                # the specified forcefield coefficients
                for atom_index1, atom_type1 in enumerate(
                    [coeff[0] for coeff in self.dpd_coeffs]
                ):
                    for atom_index2, atom_type2 in enumerate(
                        [coeff[0] for coeff in self.dpd_coeffs]
                    ):
                        self.pair_class.pair_coeff.set(
                            atom_type1,
                            atom_type2,
                            A=np.sqrt(
                                (self.dpd_coeffs[atom_index1][1] * self.e_scale)
                                * (self.dpd_coeffs[atom_index2][1] * self.e_scale)
                            ),
                            r_cut=np.sqrt(
                                (self.dpd_coeffs[atom_index1][2] * self.s_scale)
                                * (self.dpd_coeffs[atom_index2][1] * self.s_scale)
                            ),
                            gamma=self.pair_dpd_gamma_val,
                        )
                        try:
                            all_pair_types.remove(
                                "{0}-{1}".format(atom_type1, atom_type2)
                            )
                        except:
                            pass
                # Because we've been removing each pair from allPairTypes, all
                # that are left are the pair potentials that are unspecified in
                # the parXX.py (e.g. ghost particle interactions), so set these
                # interactions to zero
                for pair_type in all_pair_types:
                    self.pair_class.pair_coeff.set(
                        pair_type.split("-")[0],
                        pair_type.split("-")[1],
                        A=0.0,
                        r_cut=0.0,
                        gamma=0.0,
                    )
            elif self.pair_type.lower() == "lj":
                self.pair_class = pair.lj(r_cut=self.pair_r_cut * self.s_scale)
                self.pair_class.set_params(mode="xplor")
                for atom_index1, atom_type1 in enumerate(
                    [coeff[0] for coeff in self.lj_coeffs]
                ):
                    for atom_index2, atom_type2 in enumerate(
                        [coeff[0] for coeff in self.lj_coeffs]
                    ):
                        self.pair_class.pair_coeff.set(
                            atom_type1,
                            atom_type2,
                            epsilon=np.sqrt(
                                (self.lj_coeffs[atom_index1][1] * self.e_scale)
                                * (self.lj_coeffs[atom_index2][1] * self.e_scale)
                            ),
                            sigma=np.sqrt(
                                (self.lj_coeffs[atom_index1][2] * self.s_scale)
                                * (self.lj_coeffs[atom_index2][2] * self.s_scale)
                            ),
                        )
                        try:
                            all_pair_types.remove(
                                "{0}-{1}".format(atom_type1, atom_type2)
                            )
                        except:
                            pass
                # Because we've been removing each pair from allPairTypes, all
                # that are left are the pair potentials that are unspecified in
                # the parXX.py (e.g. ghost particle interactions), so set these
                # interactions to zero
                for pair_type in all_pair_types:
                    self.pair_class.pair_coeff.set(
                        pair_type.split("-")[0],
                        pair_type.split("-")[1],
                        epsilon=0.0,
                        sigma=0.0,
                    )
            else:
                raise SystemError(
                    "Non-dpd/lj pair potentials not yet hard-coded!"
                    " Please describe how to interpret them on this line."
                )
        # Set Bond Coeffs
        # Real bonds
        if self.bond_type.lower() == "harmonic":
            if len(self.bond_coeffs) > 0:
                self.log_quantities.append(
                    "".join(["bond_", self.bond_type, "_energy"])
                )
            self.bond_class = bond.harmonic()
            for bond_coeff in self.bond_coeffs:
                # [k] = kcal mol^{-1} \AA^{-2} * episilon/sigma^{2}, [r0] =
                # \AA * sigma^{2}
                self.bond_class.bond_coeff.set(
                    bond_coeff[0],
                    k=bond_coeff[1] * (self.e_scale / (self.s_scale ** 2)),
                    r0=bond_coeff[2] * self.s_scale,
                )
            # Ghost bonds
            # If there is no anchoring, rather than change the xml, just set the
            # bond k values to 0.
            if self.group_anchoring.lower() == "all":
                group_anchoring_types = [
                    "".join(["X", CG_type])
                    for CG_type in list(self.CG_to_template_AAIDs.keys())
                ]
            elif self.group_anchoring.lower() == "none":
                group_anchoring_types = []
            else:
                group_anchoring_types = [
                    "".join(["X", CG_type])
                    for CG_type in self.group_anchoring.split(",")
                ]
            anchor_bond_types = []
            no_anchor_bond_types = []
            for bond_type in self.AA_morphology_dict["bond"]:
                if "X" in bond_type[0]:
                    atom_type1 = bond_type[0].split("-")[0]
                    atom_type2 = bond_type[0].split("-")[1]
                    if (atom_type1 in group_anchoring_types) or (
                        atom_type2 in group_anchoring_types
                    ):
                        if bond_type[0] not in anchor_bond_types:
                            anchor_bond_types.append(bond_type[0])
                    else:
                        if bond_type[0] not in no_anchor_bond_types:
                            no_anchor_bond_types.append(bond_type[0])
            for bond_type in anchor_bond_types:
                self.bond_class.bond_coeff.set(bond_type, k=1E6, r0=0)
            for bond_type in no_anchor_bond_types:
                self.bond_class.bond_coeff.set(bond_type, k=0, r0=0)
        else:
            raise SystemError(
                "Non-harmonic bond potentials not yet hard-coded!"
                " Please describe how to interpret them on this line."
            )
        # Set Angle Coeffs
        self.angle_class = None
        if len(self.angle_coeffs) > 0:
            self.log_quantities.append("".join(["angle_", self.angle_type, "_energy"]))
            if self.angle_type.lower() == "harmonic":
                self.angle_class = angle.harmonic()
                for angle_coeff in self.angle_coeffs:
                    # [k] = kcal mol^{-1} rad^{-2} * epsilon, [t] = rad
                    self.angle_class.set_coeff(
                        angle_coeff[0],
                        k=angle_coeff[1] * self.e_scale,
                        t0=angle_coeff[2],
                    )
            else:
                raise SystemError(
                    "Non-harmonic angle potentials not yet hard-coded!"
                    " please describe how to interpret them on this line."
                )
        else:
            print("No angles detected!")
        # Set Dihedral Coeffs
        self.dihedral_class = None
        if len(self.dihedral_coeffs) > 0:
            self.log_quantities.append(
                "".join(["dihedral_", self.dihedral_type, "_energy"])
            )
            if self.dihedral_type.lower() == "table":
                self.dihedral_class = dihedral.table(width=1000)
                for dihedral_coeff in self.dihedral_coeffs:
                    self.dihedral_class.dihedral_coeff.set(
                        dihedral_coeff[0],
                        func=multi_harmonic_torsion,
                        coeff=dict(
                            v0=dihedral_coeff[1] * self.e_scale,
                            v1=dihedral_coeff[2] * self.e_scale,
                            v2=dihedral_coeff[3] * self.e_scale,
                            v3=dihedral_coeff[4] * self.e_scale,
                            v4=dihedral_coeff[5] * self.e_scale,
                        ),
                    )
            elif self.dihedral_type.lower() == "opls":
                self.dihedral_class = dihedral.opls()
                for dihedral_coeff in self.dihedral_coeffs:
                    self.dihedral_class.set_coeff(
                        dihedral_coeff[0],
                        k1=dihedral_coeff[1] * self.e_scale,
                        k2=dihedral_coeff[2] * self.e_scale,
                        k3=dihedral_coeff[3] * self.e_scale,
                        k4=dihedral_coeff[4] * self.e_scale,
                    )
            else:
                raise SystemError(
                    "Non-tabulated dihedral potentials not yet hard-coded!"
                    " Please describe how to interpret them on this line."
                )
        else:
            print("No dihedrals detected")
        # Set Improper Coeffs
        self.improper_class = None
        if len(self.improper_coeffs) > 0:
            self.improper_class = improper.harmonic()
            for improper_coeff in self.improper_coeffs:
                self.improper_class.improper_coeff.set(
                    improper_coeff[0],
                    k=improper_coeff[1] * self.e_scale,
                    chi=improper_coeff[2],
                )

    def get_integration_groups(self):
        # Based on input parameter, return all non-rigid and rigid atoms to be
        # integrated over
        if self.integration_target.lower() == "all":
            integration_types = list(self.CG_to_template_AAIDs.keys())
        else:
            integration_types = self.integration_target.split(",")
        # Add in any rigid ghost particles that might need to be integrated too
        ghost_integration_types = [
            "".join(["R", type_name]) for type_name in integration_types
        ]
        atom_IDs_to_integrate = []
        for molecule in self.CG_to_AAID_master:
            for CG_site_ID in list(molecule.keys()):
                if molecule[CG_site_ID][0] in integration_types:
                    atom_IDs_to_integrate += molecule[CG_site_ID][1]
        for atom_ID, atom_type in enumerate(self.AA_morphology_dict["type"]):
            if atom_type in ghost_integration_types:
                atom_IDs_to_integrate.append(atom_ID)
        # Create the integrateGroup which contains all of the atoms to be
        # integrated.
        # The rigidGroup constains the intersection of integrateGroup and
        # group.rigid()
        # The nonRigidGroup contains the remainder of atoms in integrateGroup
        integrate_group = group.tag_list(
            name="integrate_group", tags=atom_IDs_to_integrate
        )
        rigid_group = group.intersection(
            name="rigid_group", a=group.rigid(), b=integrate_group
        )
        non_rigid_group = group.difference(
            name="non_rigid_group", a=integrate_group, b=rigid_group
        )
        return rigid_group, non_rigid_group, integration_types


def multi_harmonic_torsion(theta, v0, v1, v2, v3, v4):
    # Definition of multiharmonic dihedral equation based on 5 input parameters
    # to be used by HOOMD
    V = (
        v0
        + v1 * np.cos(theta)
        + v2 * ((np.cos(theta)) ** 2)
        + v3 * ((np.cos(theta)) ** 3)
        + v4 * ((np.cos(theta)) ** 4)
    )
    F = (
        v1 * np.sin(theta)
        + 2 * v2 * np.cos(theta) * np.sin(theta)
        + 3 * v3 * ((np.cos(theta)) ** 2) * np.sin(theta)
        + 4 * v4 * ((np.cos(theta)) ** 3) * np.sin(theta)
    )
    return (V, F)


def obtain_scale_factors(parameter_dict):
    print("Obtaining correct scaling for epsilon and sigma...")
    # The scaling factors are 1/largestSigma in the LJ coeffs, and 1
    # / largestEpsilon
    LJFFs = []
    for CG_site, directory in parameter_dict["CG_to_template_dirs"].items():
        FF_loc = os.path.join(
            directory, parameter_dict["CG_to_template_force_fields"][CG_site]
        )
        FF = hf.load_FF_xml(FF_loc)
        LJFFs += FF["lj"]
    largest_sigma = max(list(map(float, np.array(LJFFs)[:, 2])))
    largest_epsilon = max(list(map(float, np.array(LJFFs)[:, 1])))
    return 1 / float(largest_sigma), 1 / float(largest_epsilon)


def scale_morphology(initial_morphology, parameter_dict, s_scale, e_scale):
    # If sScale != 1.0, then scale the morphology and rewrite the phase0 xml
    print("Scaling morphology by sigma = {:f}...".format(1 / s_scale))
    if s_scale != 1.0:
        hf.scale(initial_morphology, s_scale)
    hf.write_morphology_xml(
        initial_morphology,
        os.path.join(
            parameter_dict["output_morph_dir"],
            os.path.splitext(parameter_dict["morphology"])[0],
            "morphology",
            "".join(["phase0_", parameter_dict["morphology"]]),
        ),
    )


def main(
    AA_morphology_dict,
    CG_morphology_dict,
    CG_to_AAID_master,
    parameter_dict,
    chromophore_list,
):
    # Get the random seed now for all the child processes
    if parameter_dict["random_seed_override"] is not None:
        np.random.seed(parameter_dict["random_seed_override"])
    # Main execution function for run_HOOMD that performs the required MD phases
    # First, scale the input morphology based on the pair potentials such that
    # the distances and energies are normalised to the strongest pair
    # interaction and the diameter of the largest atom (makes it easier on
    # HOOMDs calculations and ensures that T = 1.0 is an interesting temperature
    # threshold)
    current_files = os.listdir(
        os.path.join(
            parameter_dict["output_morph_dir"],
            os.path.splitext(parameter_dict["morphology"])[0],
            "morphology",
        )
    )
    # sScale, eScale = obtainScaleFactors(parameterDict)
    print("Under the hood eScaling and sScaling has been disabled.")
    s_scale = 1.0
    e_scale = 1.0
    # Only scale the morphology if it hasn't been already
    if (parameter_dict["overwrite_current_data"] is False) and (
        "".join(["phase0_", parameter_dict["morphology"]]) in current_files
    ):
        pass
    else:
        scale_morphology(AA_morphology_dict, parameter_dict, s_scale, e_scale)
    # Reset logfile
    try:
        os.remove(
            os.path.join(
                parameter_dict["output_morph_dir"],
                os.path.splitext(parameter_dict["morphology"])[0],
                "morphology",
                "".join(
                    [
                        "energies_",
                        os.path.splitext(parameter_dict["morphology"])[0],
                        ".log",
                    ]
                ),
            )
        )
    except OSError:
        pass
    # Perform each molecular dynamics phase as specified in the parXX.py
    for phase_no in range(parameter_dict["number_of_phases"]):
        input_file = "phase{0:d}_{1:s}".format(phase_no, parameter_dict["morphology"])
        output_file = "phase{0:d}_{1:s}".format(
            phase_no + 1, parameter_dict["morphology"]
        )
        if output_file in current_files:
            if parameter_dict["overwrite_current_data"] is False:
                print(output_file, "already exists. Skipping...")
                continue
        md_phase(
            AA_morphology_dict,
            CG_morphology_dict,
            CG_to_AAID_master,
            parameter_dict,
            phase_no,
            os.path.join(
                parameter_dict["output_morph_dir"],
                os.path.splitext(parameter_dict["morphology"])[0],
                "morphology",
                input_file,
            ),
            os.path.join(
                parameter_dict["output_morph_dir"],
                os.path.splitext(parameter_dict["morphology"])[0],
                "morphology",
                output_file,
            ),
            s_scale,
            e_scale,
        ).optimise_structure()
    final_xml_name = os.path.join(
        parameter_dict["output_morph_dir"],
        os.path.splitext(parameter_dict["morphology"])[0],
        "morphology",
        "".join(["final_", parameter_dict["morphology"]]),
    )
    if "".join(["final_", parameter_dict["morphology"]]) not in current_files:
        # Now all phases are complete, remove the ghost particles from the
        # system
        print("Removing ghost particles to create final output...")
        remove_ghost_particles(
            os.path.join(
                parameter_dict["output_morph_dir"],
                os.path.splitext(parameter_dict["morphology"])[0],
                "morphology",
                output_file,
            ),
            final_xml_name,
            sigma=s_scale,
        )
    # Finally, update the pickle file with the most recent and realistic
    # AAMorphologyDict so that we can load it again further along the pipeline
    AA_morphology_dict = hf.load_morphology_xml(final_xml_name)
    # Now that we've obtained the final fine-grained morphology, we need to fix
    # the images to prevent issues with obtaining the chromophores and running
    # them through the ZINDO/S calculations later...
    AA_morphology_dict = hf.fix_images(AA_morphology_dict)
    # ...add in the unwrapped positions...
    AA_morphology_dict = hf.add_unwrapped_positions(AA_morphology_dict)
    # ...rewrite the final morphology xml...
    hf.write_morphology_xml(AA_morphology_dict, final_xml_name)
    # ...and write the pickle file.
    hf.write_pickle(
        (
            AA_morphology_dict,
            CG_morphology_dict,
            CG_to_AAID_master,
            parameter_dict,
            chromophore_list,
        ),
        os.path.join(
            parameter_dict["output_morph_dir"],
            os.path.splitext(parameter_dict["morphology"])[0],
            "code",
            "".join([os.path.splitext(parameter_dict["morphology"])[0], ".pickle"]),
        ),
    )
    return (
        AA_morphology_dict,
        CG_morphology_dict,
        CG_to_AAID_master,
        parameter_dict,
        chromophore_list,
    )


def remove_ghost_particles(last_phase_xml, output_file_name, sigma=1.0):
    # Remove all the ghost particles from the morphology for the final output
    final_morphology = hf.load_morphology_xml(last_phase_xml)
    # Determine the atomIDs for each particle beginning with the letters 'X'
    # or 'R' - these are the ghost particles
    atom_IDs_to_remove = []
    for atom_ID, atom_type in enumerate(final_morphology["type"]):
        if (atom_type[0] == "X") or (atom_type[0] == "R"):
            # This is a ghost particle
            atom_IDs_to_remove.append(atom_ID)
    # Reverse sorting trick so that the location indices don't change as we
    # delete particles from the system
    atom_IDs_to_remove.sort(reverse=True)
    # Now delete the atoms from the morphology
    atom_attribs = ["position", "image", "type", "mass", "diameter", "body", "charge"]
    for atom_ID in atom_IDs_to_remove:
        for key in atom_attribs:
            final_morphology[key].pop(atom_ID)
    final_morphology["natoms"] -= len(atom_IDs_to_remove)
    # Delete any constraints associated with those atoms that have been removed
    atom_constraints = ["bond", "angle", "dihedral", "improper"]
    for key in atom_constraints:
        constraints_to_remove = []
        for constraint_no, constraint in enumerate(final_morphology[key]):
            for atom_ID in constraint[1:]:
                if (atom_ID in atom_IDs_to_remove) and (
                    constraint_no not in constraints_to_remove
                ):
                    constraints_to_remove.append(constraint_no)
        constraints_to_remove.sort(reverse=True)
        for constraint_no in constraints_to_remove:
            final_morphology[key].pop(constraint_no)
    # Output the final morphology
    hf.write_morphology_xml(final_morphology, output_file_name, sigma)


if __name__ == "__main__":
    try:
        pickle_file = sys.argv[1]
    except:
        print(
            "Please specify the pickle file to load to continue the pipeline from this"
            " point."
        )
    pickle_data = hf.load_pickle(pickle_file)
    AA_morphology_dict = pickle_data[0]
    CG_morphology_dict = pickle_data[1]
    CG_to_AAID_master = pickle_data[2]
    parameter_dict = pickle_data[3]
    chromophore_list = pickle_data[4]
    main(
        AA_morphology_dict,
        CG_morphology_dict,
        CG_to_AAID_master,
        parameter_dict,
        chromophore_list,
    )
