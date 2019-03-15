import os
import pickle
import signal
import sys
import traceback
import numpy as np
import time as T
from scipy.sparse import lil_matrix
from morphct.code import helper_functions as hf


elementary_charge = 1.60217657E-19  # C
k_B = 1.3806488E-23  # m^{2} kg s^{-2} K^{-1}
hbar = 1.05457173E-34  # m^{2} kg s^{-1}
log_file = None


class carrier:
    def __init__(
        self,
        chromophore_list,
        parameter_dict,
        chromo_ID,
        lifetime,
        carrier_no,
        AA_morphology_dict,
        mol_ID_dict,
    ):
        self.ID = carrier_no
        self.image = [0, 0, 0]
        self.initial_chromophore = chromophore_list[chromo_ID]
        self.current_chromophore = chromophore_list[chromo_ID]
        if parameter_dict["hop_limit"] == 0:
            self.hop_limit = None
        else:
            self.hop_limit = parameter_dict["hop_limit"]
        self.T = parameter_dict["system_temperature"]
        self.lifetime = lifetime
        self.current_time = 0.0
        self.hole_history_matrix = None
        self.electron_history_matrix = None
        self.lambda_ij = self.current_chromophore.reorganisation_energy
        if self.current_chromophore.species.lower() == "donor":
            self.carrier_type = "hole"
            if parameter_dict["record_carrier_history"] is True:
                self.hole_history_matrix = lil_matrix(
                    (len(chromophore_list), len(chromophore_list)), dtype=int
                )
        elif self.current_chromophore.species.lower() == "acceptor":
            self.carrier_type = "electron"
            if parameter_dict["record_carrier_history"] is True:
                self.electron_history_matrix = lil_matrix(
                    (len(chromophore_list), len(chromophore_list)), dtype=int
                )
        self.no_hops = 0
        self.sim_dims = [
            [-AA_morphology_dict["lx"] / 2.0, AA_morphology_dict["lx"] / 2.0],
            [-AA_morphology_dict["ly"] / 2.0, AA_morphology_dict["ly"] / 2.0],
            [-AA_morphology_dict["lz"] / 2.0, AA_morphology_dict["lz"] / 2.0],
        ]
        self.displacement = None
        self.mol_ID_dict = mol_ID_dict
        # Set the use of average hop rates to false if the key does not exist in
        # the parameter dict
        try:
            self.use_average_hop_rates = parameter_dict["use_average_hop_rates"]
            self.average_intra_hop_rate = parameter_dict["average_intra_hop_rate"]
            self.average_inter_hop_rate = parameter_dict["average_inter_hop_rate"]
        except KeyError:
            self.use_average_hop_rates = False
        # Set the use of Koopmans' approximation to false if the key does not
        # exist in the parameter dict
        try:
            self.use_koopmans_approximation = parameter_dict[
                "use_koopmans_approximation"
            ]
        except KeyError:
            self.use_koopmans_approximation = False
        # Are we using a simple Boltzmann penalty?
        try:
            self.use_simple_energetic_penalty = parameter_dict[
                "use_simple_energetic_penalty"
            ]
        except KeyError:
            self.use_simple_energetic_penalty = False
        # Are we applying a distance penalty beyond the transfer integral?
        try:
            self.use_VRH = parameter_dict["use_VRH"]
        except KeyError:
            self.use_VRH = False
        if self.use_VRH is True:
            self.VRH_delocalisation = self.current_chromophore.VRH_delocalisation
        try:
            self.hopping_prefactor = parameter_dict["hopping_prefactor"]
        except KeyError:
            self.hopping_prefactor = 1.0

    def calculate_hop(self, chromophore_list):
        # Terminate if the next hop would be more than the termination limit
        if self.hop_limit is not None:
            if self.no_hops + 1 > self.hop_limit:
                return 1
        # Determine the hop times to all possible neighbours
        hop_times = []
        if self.use_average_hop_rates is True:
            # Use the average hop values given in the parameter dict to pick a
            # hop
            for neighbour_details in self.current_chromophore.neighbours:
                neighbour = chromophore_list[neighbour_details[0]]
                assert neighbour.ID == neighbour_details[0]
                if (
                    self.mol_ID_dict[self.current_chromophore.ID]
                    == self.mol_ID_dict[neighbour.ID]
                ):
                    hop_rate = self.average_intra_hop_rate
                else:
                    hop_rate = self.average_inter_hop_rate
                hop_time = hf.determine_event_tau(hop_rate)
                # Keep track of the chromophoreID and the corresponding tau
                hop_times.append([neighbour.ID, hop_time])
        else:
            # Obtain the reorganisation energy in J (from eV in the parameter
            # file)
            for neighbour_index, transfer_integral in enumerate(
                self.current_chromophore.neighbours_TI
            ):
                # Ignore any hops with a NoneType transfer integral (usually
                # due to an orca error)
                if transfer_integral is None:
                    continue
                delta_E_ij = self.current_chromophore.neighbours_delta_E[
                    neighbour_index
                ]
                # Load the specified hopping prefactor
                prefactor = self.hopping_prefactor
                # Get the relative image so we can update the carrier image
                # after the hop
                relative_image = self.current_chromophore.neighbours[neighbour_index][1]
                # All of the energies are in eV currently, so convert them to J
                if self.use_VRH is True:
                    neighbour_chromo = chromophore_list[
                        self.current_chromophore.neighbours[neighbour_index][0]
                    ]
                    neighbour_chromo_posn = neighbour_chromo.posn + (
                        np.array(relative_image)
                        * np.array([axis[1] - axis[0] for axis in self.sim_dims])
                    )
                    # Chromophore separation needs converting to m
                    chromophore_separation = (
                        hf.calculate_separation(
                            self.current_chromophore.posn, neighbour_chromo_posn
                        )
                        * 1E-10
                    )
                    hop_rate = hf.calculate_carrier_hop_rate(
                        self.lambda_ij * elementary_charge,
                        transfer_integral * elementary_charge,
                        delta_E_ij * elementary_charge,
                        prefactor,
                        self.T,
                        use_VRH=True,
                        rij=chromophore_separation,
                        VRH_delocalisation=self.VRH_delocalisation,
                        boltz_pen=self.use_simple_energetic_penalty,
                    )
                else:
                    hop_rate = hf.calculate_carrier_hop_rate(
                        self.lambda_ij * elementary_charge,
                        transfer_integral * elementary_charge,
                        delta_E_ij * elementary_charge,
                        prefactor,
                        self.T,
                        boltz_pen=self.use_simple_energetic_penalty,
                    )
                hop_time = hf.determine_event_tau(hop_rate)
                # Keep track of the chromophoreID and the corresponding tau
                hop_times.append(
                    [
                        self.current_chromophore.neighbours[neighbour_index][0],
                        hop_time,
                        relative_image,
                    ]
                )
        # Sort by ascending hop time
        hop_times.sort(key=lambda x: x[1])
        if len(hop_times) == 0:
            # We are trapped here, so create a dummy hop with time 1E99
            hop_times = [[self.current_chromophore.ID, 1E99, [0, 0, 0]]]
        # As long as we're not limiting by the number of hops:
        if self.hop_limit is None:
            # Ensure that the next hop does not put the carrier over its
            # lifetime
            if (self.current_time + hop_times[0][1]) > self.lifetime:
                # Send the termination signal to singleCoreRunKMC.py
                return 1
        # Move the carrier and send the contiuation signal to
        # singleCoreRunKMC.py
        # Take the quickest hop
        self.perform_hop(
            chromophore_list[hop_times[0][0]], hop_times[0][1], hop_times[0][2]
        )
        return 0

    def perform_hop(self, destination_chromophore, hop_time, relative_image):
        initial_ID = self.current_chromophore.ID
        destination_ID = destination_chromophore.ID
        self.image = list(np.array(self.image) + np.array(relative_image))
        # OLD WAY TO CALCULATE SELF.IMAGE #
        # initial_position = self.current_chromophore.posn
        # destination_position = destination_chromophore.posn
        # delta_position = destination_position - initial_position
        # for axis in range(3):
        #     half_box_length = (self.sim_dims[axis][1] - self.sim_dims[axis][0]) / 2.0
        #     while delta_position[axis] > half_box_length:
        #         # Crossed over a negative boundary, decrement image by 1
        #         delta_position[axis] -= half_box_length * 2.0
        #         self.image[axis] -= 1
        #     while delta_position[axis] < - half_box_length:
        #         # Crossed over a positive boundary, increment image by 1
        #         delta_position[axis] += half_box_length * 2.0
        #         self.image[axis] += 1
        # Carrier image now sorted, so update its current position
        self.current_chromophore = destination_chromophore
        # Increment the simulation time
        self.current_time += hop_time
        # Increment the hop counter
        self.no_hops += 1
        # Now update the sparse history matrix
        if (self.carrier_type.lower() == "hole") and (
            self.hole_history_matrix is not None
        ):
            self.hole_history_matrix[initial_ID, destination_ID] += 1
        elif (self.carrier_type.lower() == "electron") and (
            self.electron_history_matrix is not None
        ):
            self.electron_history_matrix[initial_ID, destination_ID] += 1


class termination_signal:
    kill_sent = False

    def __init__(self):
        signal.signal(signal.SIGINT, self.catch_kill)
        signal.signal(signal.SIGTERM, self.catch_kill)

    def catch_kill(self, signum, frame):
        self.kill_sent = True


class terminate(Exception):
    """This class is raised to terminate a KMC simulation"""

    def __init__(self, string):
        self.string = string

    def __str__(self):
        return self.string


def save_pickle(save_data, save_pickle_name):
    with open(save_pickle_name, "wb+") as pickle_file:
        pickle.dump(save_data, pickle_file)
    hf.write_to_file(
        log_file, ["".join(["Pickle file saved successfully as", save_pickle_name])]
    )


def calculate_displacement(initial_position, final_position, final_image, sim_dims):
    displacement = [0.0, 0.0, 0.0]
    for axis in range(3):
        displacement[axis] = (final_position[axis] - initial_position[axis]) + (
            final_image[axis] * (sim_dims[axis][1] - sim_dims[axis][0])
        )
    return np.linalg.norm(np.array(displacement))


def initialise_save_data(n_chromos, seed):
    return {
        "seed": seed,
        "ID": [],
        "image": [],
        "lifetime": [],
        "current_time": [],
        "no_hops": [],
        "displacement": [],
        "hole_history_matrix": lil_matrix((n_chromos, n_chromos), dtype=int),
        "electron_history_matrix": lil_matrix((n_chromos, n_chromos), dtype=int),
        "initial_position": [],
        "final_position": [],
        "carrier_type": [],
    }


def split_molecules(input_dictionary, chromophore_list):
    # Split the full morphology into individual molecules
    # Create a lookup table `neighbour list' for all connected atoms called
    # {bondedAtoms}
    bonded_atoms = hf.obtain_bonded_list(input_dictionary["bond"])
    molecule_list = [i for i in range(len(input_dictionary["type"]))]
    # Recursively add all atoms in the neighbour list to this molecule
    for mol_ID in range(len(molecule_list)):
        molecule_list = update_molecule(mol_ID, molecule_list, bonded_atoms)
    # Here we have a list of len(atoms) where each index gives the molID
    mol_ID_dict = {}
    for chromo in chromophore_list:
        AAID_to_check = chromo.AAIDs[0]
        mol_ID_dict[chromo.ID] = molecule_list[AAID_to_check]
    return mol_ID_dict


def update_molecule(atom_ID, molecule_list, bonded_atoms):
    # Recursively add all neighbours of atom number atomID to this molecule
    try:
        for bonded_atom in bonded_atoms[atom_ID]:
            # If the moleculeID of the bonded atom is larger than that of the
            # current one, update the bonded atom's ID to the current one's to
            # put it in this molecule, then iterate through all of the bonded
            # atom's neighbours
            if molecule_list[bonded_atom] > molecule_list[atom_ID]:
                molecule_list[bonded_atom] = molecule_list[atom_ID]
                molecule_list = update_molecule(
                    bonded_atom, molecule_list, bonded_atoms
                )
            # If the moleculeID of the current atom is larger than that of the
            # bonded one, update the current atom's ID to the bonded one's to
            # put it in this molecule, then iterate through all of the current
            # atom's neighbours
            elif molecule_list[bonded_atom] < molecule_list[atom_ID]:
                molecule_list[atom_ID] = molecule_list[bonded_atom]
                molecule_list = update_molecule(atom_ID, molecule_list, bonded_atoms)
            # Else: both the current and the bonded atom are already known to
            # be in this molecule, so we don't have to do anything else.
    except KeyError:
        # This means that there are no bonded CG sites (i.e. it's a single molecule)
        pass
    return molecule_list


def main():
    global log_file

    KMC_directory = sys.argv[1]
    CPU_rank = int(sys.argv[2])
    np.random.seed(int(sys.argv[3]))
    overwrite = False
    try:
        overwrite = bool(sys.argv[4])
    except:
        pass
    # Load `jobs_to_run' which is a list, where each element contains the
    # [carrier.ID, carrier.lifetime, carrier.carrierType]
    pickle_file_name = os.path.join(
        KMC_directory, "KMC_data_{:02d}.pickle".format(CPU_rank)
    )
    with open(pickle_file_name, "rb") as pickle_file:
        jobs_to_run = pickle.load(pickle_file)
    log_file = os.path.join(KMC_directory, "KMC_log_{:02d}.log".format(CPU_rank))
    # Reset the log file
    with open(log_file, "wb+") as log_file_handle:
        pass
    hf.write_to_file(log_file, ["Found {:d} jobs to run".format(len(jobs_to_run))])
    # Set the affinities for this current process to make sure it's maximising
    # available CPU usage
    current_PID = os.getpid()
    # try:
    #     affinity_job = sp.Popen(['taskset', '-pc', str(CPU_rank), str(current_PID)],
    #                             stdin=sp.PIPE, stdout=sp.PIPE,
    #                             stderr=sp.PIPE).communicate()
    #     # hf.write_to_file(log_file, affinity_job[0].split('\n'))
    #     # hf.write_to_file(log_file, affinity_job[1].split('\n'))
    # except OSError:
    #     hf.write_to_file(log_file, ["Taskset command not found, skipping setting of"
    #                                 " processor affinity..."])
    # Now load the main morphology pickle (used as a workaround to obtain the
    # chromophore_list without having to save it in each carrier [very memory
    # inefficient!])
    pickle_dir = KMC_directory.replace("/KMC", "/code")
    for file_name in os.listdir(pickle_dir):
        if "pickle" in file_name:
            main_morphology_pickle_name = os.path.join(pickle_dir, file_name)
    hf.write_to_file(
        log_file,
        [
            "".join(
                [
                    "Found main morphology pickle file at ",
                    main_morphology_pickle_name,
                    "! loading data...",
                ]
            )
        ],
    )
    pickle_data = hf.load_pickle(main_morphology_pickle_name)
    AA_morphology_dict = pickle_data[0]
    CG_morphology_dict = pickle_data[1]
    CG_to_AAID_master = pickle_data[2]
    parameter_dict = pickle_data[3]
    chromophore_list = pickle_data[4]
    hf.write_to_file(log_file, ["Main morphology pickle loaded!"])
    try:
        if parameter_dict["use_average_hop_rates"] is True:
            # Chosen to split hopping by inter-intra molecular hops, so get
            # molecule data
            mol_ID_dict = split_molecules(AA_morphology_dict, chromophore_list)
            # molIDDict is a dictionary where the keys are the chromoIDs, and
            # the vals are the molIDs
        else:
            raise KeyError
    except KeyError:
        mol_ID_dict = None
    # Attempt to catch a kill signal to ensure that we save the pickle before
    # termination
    killer = termination_signal()
    # Save the pickle as a list of `saveCarrier' instances that contain the
    # bare minimum
    save_data = initialise_save_data(len(chromophore_list), int(sys.argv[3]))
    if parameter_dict["record_carrier_history"] is False:
        save_data["hole_history_matrix"] = None
        save_data["electron_history_matrix"] = None
    t0 = T.time()
    save_time = T.time()
    save_slot = "slot1"
    try:
        for job_number, [carrier_no, lifetime, carrier_type] in enumerate(jobs_to_run):
            t1 = T.time()
            # Find a random position to start the carrier in
            while True:
                start_chromo_ID = np.random.randint(0, len(chromophore_list) - 1)
                if (carrier_type.lower() == "electron") and (
                    chromophore_list[start_chromo_ID].species.lower() != "acceptor"
                ):
                    continue
                elif (carrier_type.lower() == "hole") and (
                    chromophore_list[start_chromo_ID].species.lower() != "donor"
                ):
                    continue
                break
            # Create the carrier instance
            this_carrier = carrier(
                chromophore_list,
                parameter_dict,
                start_chromo_ID,
                lifetime,
                carrier_no,
                AA_morphology_dict,
                mol_ID_dict,
            )
            terminate_simulation = False
            while terminate_simulation is False:
                terminate_simulation = bool(
                    this_carrier.calculate_hop(chromophore_list)
                )
                if killer.kill_sent is True:
                    raise terminate("Kill command sent, terminating KMC simulation...")
            # Now the carrier has finished hopping, let's calculate its vitals
            initial_position = this_carrier.initial_chromophore.posn
            final_position = this_carrier.current_chromophore.posn
            final_image = this_carrier.image
            sim_dims = this_carrier.sim_dims
            this_carrier.displacement = calculate_displacement(
                initial_position, final_position, final_image, sim_dims
            )
            # Now the calculations are completed, create a barebones class
            # containing the save data
            importantData = [
                "ID",
                "image",
                "lifetime",
                "current_time",
                "no_hops",
                "displacement",
                "carrier_type",
            ]
            for name in importantData:
                save_data[name].append(this_carrier.__dict__[name])
            # Update the carrierHistoryMatrix
            if parameter_dict["record_carrier_history"] is True:
                if this_carrier.carrier_type.lower() == "hole":
                    save_data["hole_history_matrix"] += this_carrier.hole_history_matrix
                elif this_carrier.carrier_type.lower() == "electron":
                    save_data[
                        "electron_history_matrix"
                    ] += this_carrier.electron_history_matrix
            # Then add in the initial and final positions
            save_data["initial_position"].append(initial_position)
            save_data["final_position"].append(final_position)
            t2 = T.time()
            elapsed_time = float(t2) - float(t1)
            if elapsed_time < 60:
                time_units = "seconds."
            elif elapsed_time < 3600:
                elapsed_time /= 60.0
                time_units = "minutes."
            elif elapsed_time < 86400:
                elapsed_time /= 3600.0
                time_units = "hours."
            else:
                elapsed_time /= 86400.0
                time_units = "days."
            hf.write_to_file(
                log_file,
                [
                    "".join(
                        [
                            "{0:s} hopped {1:d} times, over {2:.2e} seconds, into image ".format(
                                this_carrier.carrier_type.capitalize(),
                                this_carrier.no_hops,
                                this_carrier.current_time,
                            ),
                            repr(this_carrier.image),
                            ", for a displacement of {0:.2f}, in {1:.2f} wall-clock {2:s}".format(
                                this_carrier.displacement, elapsed_time, time_units
                            ),
                        ]
                    )
                ],
            )
            # Save the pickle file every hour
            if (t2 - save_time) > 3600:
                print(
                    "Completed {0:d} of {1:d} jobs. Making checkpoint at {2:3d}%".format(
                        job_number,
                        len(jobs_to_run),
                        int(np.round((job_number + 1) / float(len(jobs_to_run)) * 100)),
                    )
                )
                hf.write_to_file(
                    log_file,
                    [
                        "Completed {0:d} of {1:d} jobs. Making checkpoint at {2:3d}%".format(
                            job_number,
                            len(jobs_to_run),
                            int(np.round((job_number + 1) / float(len(jobs_to_run)) * 100)),
                        )
                    ],
                )
                save_pickle(
                    save_data,
                    pickle_file_name.replace("data", "".join([save_slot, "_results"])),
                )
                if save_slot.lower() == "slot1":
                    save_slot = "slot2"
                elif save_slot.lower() == "slot2":
                    save_slot = "slot1"
                save_time = T.time()
    except Exception as error_message:
        print(traceback.format_exc())
        print("Saving the pickle file cleanly before termination...")
        hf.write_to_file(log_file, [str(error_message)])
        hf.write_to_file(
            log_file, ["Saving the pickle file cleanly before termination..."]
        )
        save_pickle(save_data, pickle_file_name.replace("data", "terminated_results"))
        print("Pickle saved! Exiting Python...")
        exit()
    t3 = T.time()
    elapsed_time = float(t3) - float(t0)
    if elapsed_time < 60:
        time_units = "seconds."
    elif elapsed_time < 3600:
        elapsed_time /= 60.0
        time_units = "minutes."
    elif elapsed_time < 86400:
        elapsed_time /= 3600.0
        time_units = "hours."
    else:
        elapsed_time /= 86400.0
        time_units = "days."
    hf.write_to_file(
        log_file,
        ["All jobs completed in {0:.2f} {1:s}".format(elapsed_time, time_units)],
    )
    hf.write_to_file(log_file, ["Saving the pickle file cleanly before termination..."])
    save_pickle(save_data, pickle_file_name.replace("data", "results"))
    hf.write_to_file(log_file, ["Exiting normally..."])


if __name__ == "__main__":
    main()
