import glob
import os
import pickle
import sys
import numpy as np
import subprocess as sp
from morphct.definitions import SINGLE_ORCA_RUN_FILE
from morphct.code import helper_functions as hf


class orcaError(Exception):
    def __init__(self, file_name):
        self.string = "".join(["No molecular orbital data present for ", file_name])

    def __str__(self):
        return self.string


def load_orca_output(file_name):
    with open(file_name, "r") as orca_file:
        data_file = orca_file.readlines()
    record_MO_data = False
    orbital_data = []
    for line in data_file:
        if "ORBITAL ENERGIES" in line:
            # Next line begins the MO data
            record_MO_data = True
            continue
        if record_MO_data is True:
            if "MOLECULAR ORBITALS" in line:
                # Don't need anything else from the output file
                break
            data_in_line = []
            for element in line.split(" "):
                if len(element) > 1:
                    try:
                        data_in_line.append(float(element))
                    except ValueError:
                        continue
            if len(data_in_line) == 4:
                orbital_data.append(data_in_line)
    for i in range(len(orbital_data)):
        if orbital_data[i][1] == 0:
            # This line is the first unoccupied orbital - i.e. LUMO
            LUMO = orbital_data[i][3]
            HOMO = orbital_data[i - 1][3]
            HOMO_1 = orbital_data[i - 2][3]
            LUMO_1 = orbital_data[i + 1][3]
            # Don't need any other orbitals
            break
    if record_MO_data is False:
        # Molecular orbital data not present in this file
        raise orcaError(file_name)
    return [HOMO_1, HOMO, LUMO, LUMO_1]


def modify_orca_files(file_name, failed_file, failed_count, chromophore_list):
    if failed_count == 3:
        # Three lots of reruns without any successes, try to turn off SOSCF
        print(
            "".join(
                [
                    file_name,
                    ": Three lots of reruns without any success -"
                    " turning off SOSCF to see if that helps...",
                ]
            )
        )
        turn_off_soscf(failed_file)
    elif failed_count == 6:
        # Still no joy - increase the number of SCF iterations and see if
        # convergence was just slow
        print(
            "".join(
                [
                    file_name,
                    ": Six lots of reruns without any success -"
                    " increasing the number of SCF iterations to 500...",
                ]
            )
        )
        increase_iterations(failed_file)
    elif failed_count == 9:
        # Finally, turn down the SCF tolerance
        print(
            "".join(
                [
                    file_name,
                    ": Nine lots of reruns without any success -"
                    " decreasing SCF tolerance (sloppySCF)...",
                ]
            )
        )
        reduce_tolerance(failed_file)
    elif failed_count == 12:
        print(
            "".join(
                [
                    file_name,
                    ": Failed to rerun orca 12 times, one final thing"
                    " that can be done is to change the numerical accuracy...",
                ]
            )
        )
        revert_orca_files(failed_file)
        increase_grid(failed_file)
    elif failed_count == 15:
        print(
            "".join(
                [
                    file_name,
                    ": Failed to rerun orca 15 times. Will try high"
                    " numerical accuracy with no SOSCF as a last-ditch effort...",
                ]
            )
        )
        increase_grid_no_soscf(failed_file)
    elif failed_count == 18:
        # SERIOUS PROBLEM
        print(
            "".join(
                [
                    file_name,
                    ": Failed to rerun orca 18 times, even with all"
                    " the input file tweaks. Examine the geometry - it is most likely"
                    " unreasonable.",
                ]
            )
        )
        file_string = os.path.splitext(os.path.split(file_name)[1])[0]
        for chromo_string in file_string.split("-"):
            chromo_ID = int(chromo_string)
            print("AAIDs for chromophore", chromo_ID)
            print(chromophore_list[chromo_ID].AAIDs)
        print("Reverting {:s} back to its original state...".format(file_name))
        revert_orca_files(failed_file)
        return 1
    return 0


def turn_off_soscf(input_file):
    with open(input_file, "r") as file_name:
        original_lines = file_name.readlines()
    original_lines[3] = "!ZINDO/S NOSOSCF\n"
    with open(input_file, "w+") as file_name:
        file_name.writelines(original_lines)


def reduce_tolerance(input_file):
    with open(input_file, "r") as file_name:
        original_lines = file_name.readlines()
    original_lines[3] = "!ZINDO/S NoSOSCF SloppySCF\n"
    with open(input_file, "w+") as file_name:
        file_name.writelines(original_lines)


def increase_iterations(input_file):
    with open(input_file, "r") as file_name:
        original_lines = file_name.readlines()
    original_lines.append("\n%scf MaxIter 500 end")
    with open(input_file, "w+") as file_name:
        file_name.writelines(original_lines)


def increase_grid(input_file):
    with open(input_file, "r") as file_name:
        original_lines = file_name.readlines()
    original_lines[3] = "!ZINDO/S SlowConv Grid7 NoFinalGrid\n"
    with open(input_file, "w+") as file_name:
        file_name.writelines(original_lines)


def increase_grid_no_soscf(input_file):
    with open(input_file, "r") as file_name:
        original_lines = file_name.readlines()
    original_lines[3] = "!ZINDO/S SlowConv Grid7 NoFinalGrid NoSOSCF SloppySCF\n"
    original_lines.append("\n%scf MaxIter 500 end")
    with open(input_file, "w+") as file_name:
        file_name.writelines(original_lines)


def revert_orca_files(input_file):
    with open(input_file, "r") as file_name:
        original_lines = file_name.readlines()
    original_lines[3] = "! ZINDO/S\n"
    for line_no in range(len(original_lines)):
        # REMOVE THE SCF ITER
        if "%scf MaxIter" in original_lines[line_no]:
            original_lines.pop(line_no)
            break
    with open(input_file, "w+") as file_name:
        file_name.writelines(original_lines)


def rerun_fails(failed_chromo_files, parameter_dict, chromophore_list):
    print("")
    print(failed_chromo_files)
    number_of_fails = len(list(failed_chromo_files.keys()))
    if number_of_fails == 1:
        print("There was 1 failed job.")
    else:
        print("There were {:d} failed jobs.".format(number_of_fails))
    proc_IDs = parameter_dict["proc_IDs"]
    pop_list = []
    permanently_failed = {}
    # Firstly, modify the input files to see if numerical tweaks make orca
    # happier
    for failed_file, failed_data in failed_chromo_files.items():
        failed_count = failed_data[0]
        error_code = modify_orca_files(
            failed_file,
            os.path.join(
                parameter_dict["output_orca_directory"],
                "chromophores",
                "input_orca",
                failed_file.replace(".out", ".inp"),
            ),
            failed_count,
            chromophore_list,
        )
        if error_code == 1:
            # Don't delete the elements from the list here because we're still
            # trying to iterate over this dict and it cannot change length!
            pop_list.append(failed_file)
            permanently_failed[failed_file] = failed_data
    # Now pop the correct elements from the failed_chromo_files dict
    for failed_file in pop_list:
        failed_chromo_files.pop(failed_file)
    # If there are no files left, then everything has failed so this function
    # has completed its task
    if len(failed_chromo_files) == 0:
        return failed_chromo_files, permanently_failed
    # Otherwise, rerun those failed files.
    # First, find the correct locations of the input Files
    input_files = [
        os.path.join(
            parameter_dict["output_orca_directory"],
            "chromophores",
            "input_orca",
            file_name.replace(".out", ".inp"),
        )
        for file_name in list(failed_chromo_files.keys())
    ]
    # As before, split the list of reruns based on the number of processors
    jobs_list = [
        input_files[i : i + (int(np.ceil(len(input_files) / len(proc_IDs)))) + 1]
        for i in range(
            0, len(input_files), int(np.ceil(len(input_files) / float(len(proc_IDs))))
        )
    ]
    print(jobs_list)
    # Write the jobs pickle for single_core_run_orca to obtain
    with open(
        os.path.join(
            parameter_dict["output_orca_directory"], "chromophores", "orca_jobs.pickle"
        ),
        "wb+",
    ) as pickle_file:
        pickle.dump(jobs_list, pickle_file)
    # Now rerun orca
    if len(jobs_list) <= len(proc_IDs):
        proc_IDs = proc_IDs[: len(jobs_list)]
    running_jobs = []
    for CPU_rank in proc_IDs:
        # The final argument here tells orca to ignore the presence of the
        # output file and recalculate
        running_jobs.append(
            sp.Popen(
                [
                    "python",
                    SINGLE_ORCA_RUN_FILE,
                    parameter_dict["output_orca_directory"],
                    parameter_dict["output_morphology_directory"],
                    str(CPU_rank),
                    "1",
                    "0",
                ]
            )
        )
    # Wait for running jobs to finish
    [p.wait() for p in running_jobs]
    # Finally, return the failed files list to the main failure handler to see
    # if we need to iterate
    return failed_chromo_files, permanently_failed


def calculate_delta_E(chromophore_list, chromo1_ID, chromo2_ID):
    chromo1 = chromophore_list[chromo1_ID]
    chromo2 = chromophore_list[chromo2_ID]
    # NOTE: SANITY CHECK
    if (chromo1.ID != chromo1_ID) or (chromo2.ID != chromo2_ID):
        raise SystemError(
            "chromo1.ID ({0:d}) != chromo1_ID ({1:d}), or chromo2.ID ({2:d}) != chromo2_ID ({3:d})! CHECK CODE!".format(
                chromo1.ID, chromo1_ID, chromo2.ID, chromo2_ID
            )
        )
    # END OF SANITY CHECK
    if chromo1.species.lower() == "donor":
        # Hole transporter
        chromo1_E = chromo1.HOMO
    elif chromo1.species.lower() == "acceptor":
        # Electron transporter
        chromo1_E = chromo1.LUMO
    if chromo2.species.lower() == "donor":
        # Hole transporter
        chromo2_E = chromo2.HOMO
    elif chromo2.species.lower() == "acceptor":
        # Electron transporter
        chromo2_E = chromo2.LUMO
    return chromo2_E - chromo1_E


def calculate_TI(orbital_splitting, delta_E):
    # Use the energy splitting in dimer method to calculate the electronic
    # transfer integral in eV
    if delta_E ** 2 > orbital_splitting ** 2:
        # Avoid an imaginary TI by returning zero.
        # (Could use KOOPMAN'S APPROXIMATION here if desired)
        TI = 0
    else:
        TI = 0.5 * np.sqrt((orbital_splitting ** 2) - (delta_E ** 2))
    return TI


def update_single_chromophore_list(chromophore_list, parameter_dict):
    orca_input_dir = os.path.join(
        parameter_dict["output_orca_directory"], "chromophores", "input_orca"
    )
    orca_output_dir = os.path.join(
        parameter_dict["output_orca_directory"], "chromophores", "output_orca"
    )
    # NOTE: This can possibly be done by recursively iterating through the
    # neighbourlist of each chromophore, but I imagine Python will whinge about
    # the levels of recursion, so for now I'll just go through every
    # chromophore twice.
    # Firstly, set the energy levels for each single chromophore, rerunning
    # them if they fail.
    # failed_single_chromos has the form {'file_name': [fail_count,
    # location_in_chromophore_list]}
    failed_single_chromos = {}
    for chromo_location, chromophore in enumerate(chromophore_list):
        file_name = "single/{:05d}.out".format(chromophore.ID)
        print("\rDetermining energy levels for", file_name, end=" ")
        if sys.stdout is not None:
            sys.stdout.flush()
        # Update the chromophores in the chromophore_list with their
        # energy_levels
        try:
            energy_levels = load_orca_output(os.path.join(orca_output_dir, file_name))
            chromophore.HOMO_1 = energy_levels[0]
            chromophore.HOMO = energy_levels[1]
            chromophore.LUMO = energy_levels[2]
            chromophore.LUMO_1 = energy_levels[3]
            # If this file had originally failed, then we can safely remove it
            # from the fail list.
            if file_name in failed_single_chromos.keys():
                failed_single_chromos.pop(file_name)
        except orcaError:
            failed_single_chromos[file_name] = [1, chromo_location]
            continue
    print("")
    # Rerun any failed orca jobs
    while len(failed_single_chromos) > 0:
        failed_single_chromos, permanently_failed = rerun_fails(
            failed_single_chromos, parameter_dict, chromophore_list
        )
        if len(permanently_failed) > 0:
            print(permanently_failed)
            print("--== CRITICAL ERROR ==--")
            print(
                "THE ABOVE SINGLE-CHROMOPHORE SYSTEMS FAILED PERMANENTLY. THESE NEED"
                " FIXING/REMOVING FROM THE SYSTEM BEFORE ANY FURTHER DATA CAN BE"
                " OBTAINED."
            )
            if parameter_dict["remove_orca_inputs"]:
                print("Deleting the remaining orca inputs...")
                shutil.rmtree(orca_input_dir)
            if parameter_dict["remove_orca_outputs"]:
                print("Deleting the remaining orca inputs...")
                shutil.rmtree(orca_output_dir)
            exit()
        successful_reruns = []
        # Now check all of the files to see if we can update the
        # chromophore_list
        for chromo_name, chromo_data in failed_single_chromos.items():
            print("Checking previously failed", chromo_name)
            chromo_ID = chromo_data[1]
            try:
                # Update the chromophore data in the chromophore_list
                energy_levels = load_orca_output(
                    os.path.join(orca_output_dir, chromo_name)
                )
                chromophore_list[chromo_ID].HOMO_1 = energy_levels[0]
                chromophore_list[chromo_ID].HOMO = energy_levels[1]
                chromophore_list[chromo_ID].LUMO = energy_levels[2]
                chromophore_list[chromo_ID].LUMO_1 = energy_levels[3]
                # This chromophore didn't fail, so remove it from the failed
                # list
                successful_reruns.append(chromo_name)
            except orcaError:
                # This chromophore failed so increment its fail counter
                failed_single_chromos[chromo_name][0] += 1
                continue
        for chromo_name in successful_reruns:
            failed_single_chromos.pop(chromo_name)
    print("")
    return chromophore_list


def update_pair_chromophore_list(chromophore_list, parameter_dict):
    # Now that all the single chromophore energy levels are done, iterate
    # through again and check the neighbours, rerunning the pair file if it
    # failed (which it won't have done because all my chromophores are
    # delicious now).
    orca_input_dir = os.path.join(
        parameter_dict["output_orca_directory"], "chromophores", "input_orca"
    )
    orca_output_dir = os.path.join(
        parameter_dict["output_orca_directory"], "chromophores", "output_orca"
    )
    failed_pair_chromos = {}
    for chromo_location, chromophore in enumerate(chromophore_list):
        neighbour_IDs = [neighbour_data[0] for neighbour_data in chromophore.neighbours]
        for neighbour_loc, neighbour_ID in enumerate(neighbour_IDs):
            if chromophore.ID > neighbour_ID:
                continue
            file_name = "pair/{0:05d}-{1:05d}.out".format(chromophore.ID, neighbour_ID)
            print("\rDetermining energy levels for", file_name, end=" ")
            if sys.stdout is not None:
                sys.stdout.flush()
            try:
                energy_levels = load_orca_output(
                    os.path.join(orca_output_dir, file_name)
                )
                dimer_HOMO_1 = energy_levels[0]
                dimer_HOMO = energy_levels[1]
                dimer_LUMO = energy_levels[2]
                dimer_LUMO_1 = energy_levels[3]
                # If this file had originally failed, then we can safely remove it
                # from the fail list.
                if file_name in failed_pair_chromos.keys():
                    failed_pair_chromos.pop(file_name)
            except orcaError:
                failed_pair_chromos[file_name] = [1, chromo_location, neighbour_ID]
                continue
            # Calculate the delta_E between the two single chromophores
            try:
                if parameter_dict["use_koopmans_approximation"]:
                    delta_E = 0.0
                else:
                    # Calculate Delta_E normally
                    raise KeyError
            except KeyError:
                delta_E = calculate_delta_E(
                    chromophore_list, chromophore.ID, neighbour_ID
                )
            # Check the chromophore species
            assert (
                chromophore_list[chromophore.ID].species
                == chromophore_list[neighbour_ID].species
            )
            species = chromophore_list[chromophore.ID].species
            # Calculate the TI using the ESD method
            if species.lower() == "donor":
                TI = calculate_TI(dimer_HOMO - dimer_HOMO_1, delta_E)
            elif species.lower() == "acceptor":
                TI = calculate_TI(dimer_LUMO - dimer_LUMO_1, delta_E)
            # Get the location of the current chromophore.ID in the neighbour's
            # neighbourList
            reverse_loc = [
                neighbour_data[0]
                for neighbour_data in chromophore_list[neighbour_ID].neighbours
            ].index(chromophore.ID)
            # Update both the current chromophore and the neighbour (for the
            # reverse hop)
            chromophore.neighbours_delta_E[neighbour_loc] = delta_E
            chromophore_list[neighbour_ID].neighbours_delta_E[reverse_loc] = -delta_E
            chromophore.neighbours_TI[neighbour_loc] = TI
            chromophore_list[neighbour_ID].neighbours_TI[reverse_loc] = TI
            # DEBUG ASSERTIONS
            # Check list index corresponds to chromophore ID
            assert chromo_location == chromophore_list[chromo_location].ID
            assert chromo_location == chromophore.ID
            # Check the neighbourLoc and reverseLoc give the correct
            # chromophoreIDs
            assert (
                chromophore_list[chromophore.ID].neighbours[neighbour_loc][0]
                == chromophore_list[neighbour_ID].ID
            )
            assert (
                chromophore_list[neighbour_ID].neighbours[reverse_loc][0]
                == chromophore_list[chromophore.ID].ID
            )
            # Check the chromophoreList has been updated after updating the
            # chromophore instance
            assert (
                chromophore_list[chromophore.ID].neighbours_TI[neighbour_loc]
                == chromophore.neighbours_TI[neighbour_loc]
            )
            # Check the TI of the forward and backward hops are the same
            assert (
                chromophore_list[chromophore.ID].neighbours_TI[neighbour_loc]
                == chromophore_list[neighbour_ID].neighbours_TI[reverse_loc]
            )
            # Check the chromophoreList has been updated after updating the
            # chromophore instance
            assert (
                chromophore_list[chromophore.ID].neighbours_delta_E[neighbour_loc]
                == chromophore.neighbours_delta_E[neighbour_loc]
            )
            # Check the Delta_E of the forward and backward hops are *= -1
            assert (
                chromophore_list[chromophore.ID].neighbours_delta_E[neighbour_loc]
                == -chromophore_list[neighbour_ID].neighbours_delta_E[reverse_loc]
            )
            # END DEBUG ASSERTIONS
    print("")
    while len(failed_pair_chromos) > 0:
        failed_pair_chromos, permanently_failed = rerun_fails(
            failed_pair_chromos, parameter_dict, chromophore_list
        )
        if len(permanently_failed) > 0:
            print("--== WARNING ==--")
            print(
                "The above chromophore-pair systems failed permanently. Setting their"
                " transfer integrals to zero, preventing these hops from ever taking"
                " place in the KMC."
            )
            for file_name, chromo_data in permanently_failed.items():
                chromo1_ID = chromo_data[1]
                chromo2_ID = chromo_data[2]
                TI = 0.0
                delta_E = 0.0
                # Get the location of the neighbour's ID in the current
                # chromophores's neighbourList
                neighbour_loc = [
                    neighbour_data[0]
                    for neighbour_data in chromophore_list[chromo1_ID].neighbours
                ].index(chromo2_ID)
                # Get the location of the current chromophore's ID in the
                # neighbour's neighbourList
                reverse_loc = [
                    neighbour_data[0]
                    for neighbour_data in chromophore_list[chromo2_ID].neighbours
                ].index(chromo1_ID)
                # Update both the current chromophore and the neighbour (for the reverse
                # hop)
                chromophore_list[chromo1_ID].neighbours_delta_E[neighbour_loc] = delta_E
                chromophore_list[chromo2_ID].neighbours_delta_E[reverse_loc] = -delta_E
                chromophore_list[chromo1_ID].neighbours_TI[neighbour_loc] = TI
                chromophore_list[chromo2_ID].neighbours_TI[reverse_loc] = TI
        successful_reruns = []
        for file_name, chromo_data in failed_pair_chromos.items():
            print("Checking previously failed", file_name)
            chromo1_ID = chromo_data[1]
            chromo2_ID = chromo_data[2]
            try:
                energy_levels = load_orca_output(
                    os.path.join(orca_output_dir, file_name)
                )
                dimer_HOMO_1 = energy_levels[0]
                dimer_HOMO = energy_levels[1]
                dimer_LUMO = energy_levels[2]
                dimer_LUMO_1 = energy_levels[3]
            except orcaError:
                # This dimer failed so increment its fail counter
                failed_pair_chromos[file_name][0] += 1
                print(file_name, "still failed, incrementing counter")
                continue
            # Calculate the delta_E between the two single chromophores
            try:
                if parameter_dict["use_koopmans_approximation"]:
                    delta_E = 0.0
                else:
                    # Calculate Delta_E normally
                    raise KeyError
            except KeyError:
                delta_E = calculate_delta_E(
                    chromophore_list, chromophore.ID, neighbour_ID
                )
            # Check the chromophore species
            assert (
                chromophore_list[chromophore.ID].species
                == chromophore_list[neighbour_ID].species
            )
            species = chromophore_list[chromophore.ID].species
            # Calculate the TI using the ESD method
            if species.lower() == "donor":
                TI = calculate_TI(dimer_HOMO - dimer_HOMO_1, delta_E)
            elif species.lower() == "acceptor":
                TI = calculate_TI(dimer_LUMO - dimer_LUMO_1, delta_E)
            # Get the location of the neighbour's ID in the current
            # chromophores's neighbourList
            neighbour_loc = [
                neighbour_data[0]
                for neighbour_data in chromophore_list[chromo1_ID].neighbours
            ].index(chromo2_ID)
            # Get the location of the current chromophore's ID in the
            # neighbour's neighbourList
            reverse_loc = [
                neighbour_data[0]
                for neighbour_data in chromophore_list[chromo2_ID].neighbours
            ].index(chromo1_ID)
            # Update both the current chromophore and the neighbour (for the
            # reverse hop)
            chromophore_list[chromo1_ID].neighbours_delta_E[neighbour_loc] = delta_E
            chromophore_list[chromo2_ID].neighbours_delta_E[reverse_loc] = -delta_E
            chromophore_list[chromo1_ID].neighbours_TI[neighbour_loc] = TI
            chromophore_list[chromo2_ID].neighbours_TI[reverse_loc] = TI
            # This rerun was successful so remove this chromophore from the
            # rerun list
            successful_reruns.append(file_name)
            print(file_name, "was successful!")
            # DEBUG ASSERTIONS
            # Check the neighbourLoc and reverseLoc give the correct
            # chromophoreIDs
            assert (
                chromophore_list[chromo1_ID].neighbours[neighbour_loc][0]
                == chromophore_list[chromo2_ID].ID
            )
            assert (
                chromophore_list[chromo2_ID].neighbours[reverse_loc][0]
                == chromophore_list[chromo1_ID].ID
            )
            # Check the TI of the forward and backward hops are the same
            assert (
                chromophore_list[chromo1_ID].neighbours_TI[neighbour_loc]
                == chromophore_list[chromo2_ID].neighbours_TI[reverse_loc]
            )
            # Check the Delta_E of the forward and backward hops are *= -1
            assert (
                chromophore_list[chromo1_ID].neighbours_delta_E[neighbour_loc]
                == -chromophore_list[chromo2_ID].neighbours_delta_E[reverse_loc]
            )
            # END DEBUG ASSERTIONS
        for file_name in successful_reruns:
            failed_pair_chromos.pop(file_name)
    print("")
    return chromophore_list


def scale_energies(chromophore_list, parameter_dict):
    # Shorter chromophores have significantly deeper HOMOs because they are
    # treated as small molecules instead of chain segments. To rectify this,
    # find the average energy level for each chromophore and then map that
    # average to the literature value.
    # First, get the energy level data
    chromophore_species = {k: [] for k in parameter_dict["chromophore_species"].keys()}
    chromophore_MO_info = {k: {} for k in parameter_dict["chromophore_species"].keys()}
    for chromo in chromophore_list:
        chromophore_species[chromo.sub_species].append(chromo.get_MO_energy())

    for sub_species, chromo_energy in chromophore_species.items():
        lit_DOS_std = parameter_dict["chromophore_species"][sub_species][
            "target_DOS_std"
        ]
        lit_MO = parameter_dict["chromophore_species"][sub_species]["literature_MO"]
        chromophore_MO_info[sub_species]["target_DOS_std"] = lit_DOS_std
        chromophore_MO_info[sub_species]["av_MO"] = np.average(chromo_energy)
        chromophore_MO_info[sub_species]["std_MO"] = np.std(chromo_energy)
        chromophore_MO_info[sub_species]["E_shift"] = (
            lit_MO - chromophore_MO_info[sub_species]["av_MO"]
        )

    for chromo in chromophore_list:
        E_shift = chromophore_MO_info[chromo.sub_species]["E_shift"]
        target_DOS_std = chromophore_MO_info[chromo.sub_species]["target_DOS_std"]
        std_MO = chromophore_MO_info[chromo.sub_species]["std_MO"]
        av_MO = chromophore_MO_info[chromo.sub_species]["av_MO"]

        chromo.HOMO_1 += E_shift
        chromo.HOMO += E_shift
        chromo.LUMO += E_shift
        chromo.LUMO_1 += E_shift

        if (target_DOS_std is not None) and (target_DOS_std < std_MO):
            # Determine how many sigmas away from the mean this datapoint is
            sigma = (chromo.get_MO_energy() - av_MO) / std_MO
            # Calculate the new deviation from the mean based on the target
            # STD and sigma
            new_deviation = target_DOS_std * sigma
            # Work out the change in energy to be applied to meet this target
            # energy level
            delta_E = (av_MO + new_deviation) - chromo.get_MO_energy()
            # Apply the energy level displacement
            chromo.HOMO_1 += delta_E
            chromo.HOMO += delta_E
            chromo.LUMO += delta_E
            chromo.LUMO_1 += delta_E
    return chromophore_list


def main(
    AA_morphology_dict,
    CG_morphology_dict,
    CG_to_AAID_master,
    parameter_dict,
    chromophore_list,
):
    pickle_name = os.path.join(
        parameter_dict["output_morphology_directory"],
        "code",
        "".join([os.path.splitext(parameter_dict["morphology"])[0], ".pickle"]),
    )
    # First, check that we need to examine the single chromophores
    run_singles = False
    if parameter_dict["overwrite_current_data"] is False:
        # Only perform this check if the user hasn't already specified to
        # overwrite the data (in which case it runs anyway)
        # Run all singles if any of the single's data is missing (i.e. the
        # HOMO level should suffice because all energy levels are updated at
        # the same time, so we don't need to check all of them individually)
        for chromophore in chromophore_list:
            if chromophore.HOMO is None:
                run_singles = True
    if (run_singles is True) or (parameter_dict["overwrite_current_data"] is True):
        print("Beginning analysis of single chromophores...")
        chromophore_list = update_single_chromophore_list(
            chromophore_list, parameter_dict
        )
        # Now include any scaling to narrow the DoS or modulate the mean to
        # match the literature HOMO/LUMO levels (which helps to negate the
        # effect of short chromophores with additional hydrogens/terminating
        # groups)
        print("Scaling energies...")
        chromophore_list = scale_energies(chromophore_list, parameter_dict)
        print("Single chromophore calculations completed. Saving...")
        hf.write_pickle(
            (
                AA_morphology_dict,
                CG_morphology_dict,
                CG_to_AAID_master,
                parameter_dict,
                chromophore_list,
            ),
            pickle_name,
        )
        print("Deleting outputs...")
        if parameter_dict["remove_orca_outputs"] is True:
            for chromophore in chromophore_list:
                try:
                    os.remove(os.path.join(orca_output_dir, chromophore.orca_output))
                except FileNotFoundError:
                    # Already deleted
                    pass
    else:
        print("All single chromophore calculations already performed. Skipping...")
    # Then, check the pairs
    run_pairs = False
    if parameter_dict["overwrite_current_data"] is False:
        for chromophore in chromophore_list:
            # Just check the first neighbour for each chromophore
            for neighbour in chromophore.neighbours_TI:
                if neighbour is None:
                    run_pairs = True
                    break
    if (run_pairs is True) or (parameter_dict["overwrite_current_data"] is True):
        print("Beginning analysis of chromophore pairs...")
        chromophore_list = update_pair_chromophore_list(
            chromophore_list, parameter_dict
        )
        # DEBUG Testing - you can remove these as the assertions in
        # update_pair_chromophore_list should already cover them, however they
        # are fast and will ensure that there are no errors in the
        # chromophore_list after calculating the T_ij and Delta_E_ijs
        T_ij_error = check_forward_backward_hop_T_ij(chromophore_list)
        delta_E_error = check_forward_backward_hop_E_ij(chromophore_list)
        if T_ij_error or delta_E_error:
            raise SystemError("assertions failed, please address in code.")
        # END OF DEBUG Testing
        print("Pair chromophore calculations completed. Saving...")
        hf.write_pickle(
            (
                AA_morphology_dict,
                CG_morphology_dict,
                CG_to_AAID_master,
                parameter_dict,
                chromophore_list,
            ),
            pickle_name,
        )
        if parameter_dict["remove_orca_outputs"] is True:
            for chromophore in chromophore_list:
                neighbour_IDs = [
                    neighbour_data[0] for neighbour_data in chromophore.neighbours
                ]
                for neighbour_loc, neighbour_ID in enumerate(neighbour_IDs):
                    file_name = "pair/{0:05d}-{1:05d}.out".format(
                        chromophore.ID, neighbour_ID
                    )
                try:
                    os.remove(os.path.join(orca_output_dir, file_name))
                except FileNotFoundError:
                    # Already deleted
                    pass
    else:
        print("All pair chromophore calculations already performed. Skipping...")
    return (
        AA_morphology_dict,
        CG_morphology_dict,
        CG_to_AAID_master,
        parameter_dict,
        chromophore_list,
    )


def check_forward_backward_hop_T_ij(chromophore_list):
    # Check reverse lookup: T_ij === T_ji
    donor_errors = 0
    acceptor_errors = 0
    for chromo1 in chromophore_list:
        chromo1_ID = chromo1.ID
        for neighbour_index, chromo2_details in enumerate(chromo1.neighbours):
            chromo2_ID = chromo2_details[0]
            chromo1_to_2_TI = chromo1.neighbours_TI[neighbour_index]
            # Sanity check
            assert chromo2_ID == chromophore_list[chromo2_ID].ID
            chromo2 = chromophore_list[chromo2_ID]
            neighbour2_index = 0
            for neighbour2_index, chromo1_details in enumerate(chromo2.neighbours):
                if chromo1_details[0] != chromo1_ID:
                    continue
                chromo2_to_1_TI = chromo2.neighbours_TI[neighbour2_index]
                break
            assert chromo1.species == chromo2.species
            try:
                assert chromo1_to_2_TI == chromo2_to_1_TI
                # Put other assertions in here
            except AssertionError:
                print("\n<ERROR FOUND>")
                print("1 to 2", chromo1_to_2_TI)
                print("2 to 1", chromo2_to_1_TI)
                print("Chromo 1 ID =", chromo1_ID, "Chromo 2 ID =", chromo2_ID)
                print(
                    "Chromo1 Neighbours: Look for index =",
                    neighbour_index,
                    "in",
                    chromo1.neighbours,
                )
                print(
                    "Chromo2 Neighbours: Look for index =",
                    neighbour2_index,
                    "in",
                    chromo2.neighbours,
                )
                print(
                    "Chromo1 TIs: Look for index =",
                    neighbour_index,
                    "in",
                    chromo1.neighboursTI,
                )
                print(
                    "Chromo2 TIs: Look for index =",
                    neighbour2_index,
                    "in",
                    chromo2.neighboursTI,
                )
                if chromo1.species.lower() == "donor":
                    donor_errors += 1
                elif chromo1.species.lower() == "acceptor":
                    acceptor_errors += 1
    if (donor_errors > 0) or (acceptor_errors > 0):
        print("--== CRITICAL ERROR ==--")
        print(
            "\nThere were",
            donor_errors,
            "cases where T_ij != T_ji in the donor chromophores.",
        )
        print(
            "\nThere were",
            acceptor_errors,
            "cases where T_ij != T_ji in the acceptor chromophores.",
        )
        return 1
    return 0


def check_forward_backward_hop_E_ij(chromophore_list):
    # Check reverse lookup: Delta E_ij === -Delta E_ji
    donor_errors = 0
    acceptor_errors = 0
    for chromophore in chromophore_list:
        chromo_ID = chromophore.ID
        for neighbour_loc, neighbour_deets in enumerate(chromophore.neighbours):
            neighbour_ID = neighbour_deets[0]
            assert chromo_ID == chromophore_list[chromo_ID].ID
            assert neighbour_ID == chromophore_list[neighbour_ID].ID
            # Get the location of the current chromophore.ID in the neighbour's
            # neighbourList
            reverse_loc = [
                neighbour_data[0]
                for neighbour_data in chromophore_list[neighbour_ID].neighbours
            ].index(chromophore.ID)
            assert (
                neighbour_ID == chromophore_list[chromo_ID].neighbours[neighbour_loc][0]
            )
            assert (
                chromo_ID == chromophore_list[neighbour_ID].neighbours[reverse_loc][0]
            )
            # Update both the current chromophore and the neighbour (for the
            # reverse hop)
            try:
                assert (
                    chromophore_list[chromo_ID].neighbours_delta_E[neighbour_loc]
                    == -chromophore_list[neighbour_ID].neighbours_delta_E[reverse_loc]
                )
            except AssertionError:
                print("\nHOP FROM", chromo_ID, "TO", neighbour_ID)
                print(
                    neighbour_ID,
                    "should be here",
                    chromophore.neighbours[neighbour_loc],
                )
                print(
                    chromo_ID,
                    "should be here",
                    chromophore_list[neighbour_ID].neighbours[reverse_loc],
                )
                print("--== Transfer Integrals ==--")
                print(
                    "FORWARD:",
                    chromophore_list[chromo_ID].neighbours_TI[neighbour_loc],
                    "backward:",
                    chromophore_list[neighbour_ID].neighbours_TI[reverse_loc],
                )
                print("--== Delta E_ij ==--")
                print(
                    "FORWARD:",
                    chromophore_list[chromo_ID].neighbours_delta_E[neighbour_loc],
                    "backward:",
                    chromophore_list[neighbour_ID].neighbours_delta_E[reverse_loc],
                )
                if chromophore.species.lower() == "donor":
                    donor_errors += 1
                elif chromophore.species.lower() == "acceptor":
                    acceptor_errors += 1
    if (donor_errors > 0) or (acceptor_errors > 0):
        print("--== CRITICAL ERROR ==--")
        print(
            "\nThere were",
            donor_errors,
            "cases where E_ij != -E_ji in the donor chromophores.",
        )
        print(
            "\nThere were",
            acceptor_errors,
            "cases where E_ij != -E_ji in the acceptor chromophores.",
        )
        return 1
    return 0


if __name__ == "__main__":
    try:
        pickle_file = sys.argv[1]
    except NameError:
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
