import os
import pickle
import sys
import multiprocessing as mp
import numpy as np
import subprocess as sp
from morphct.definitions import PROJECT_ROOT, SINGLE_ORCA_RUN_FILE
from morphct.code import helper_functions as hf


def create_input_files(chromophore_list, AA_morphology_dict, parameter_dict):
    # Singles first
    for chromophore in chromophore_list:
        # Include the molecule terminating units on the required atoms of the
        # chromophore
        if chromophore.terminate is True:
            terminating_group_positions = terminate_monomers(
                chromophore, parameter_dict, AA_morphology_dict
            )
            terminating_group_images = [chromophore.image] * len(
                terminating_group_positions
            )
        else:
            terminating_group_positions = None
            terminating_group_images = None
        write_orca_inp(
            AA_morphology_dict,
            chromophore.AAIDs,
            [chromophore.image] * len(chromophore.AAIDs),
            terminating_group_positions,
            terminating_group_images,
            "".join([parameter_dict["output_orca_directory"], chromophore.orca_input]),
        )
    print("")
    # Determine how many pairs there are first:
    number_of_pairs = 0
    for chromo in chromophore_list:
        for neighbour in chromo.neighbours:
            if int(neighbour[0]) > chromo.ID:
                number_of_pairs += 1
    number_of_pairs = np.sum([len(chromo.neighbours) for chromo in chromophore_list])
    print("There are", int(number_of_pairs / 2), "total neighbour pairs to consider.")
    # /2 because the forwards and backwards hops are identical
    # Then consider each chromophore against every other chromophore
    for chromophore1 in chromophore_list:
        neighbours_ID = [neighbour[0] for neighbour in chromophore1.neighbours]
        neighbours_image = [neighbour[1] for neighbour in chromophore1.neighbours]
        for chromophore2 in chromophore_list:
            # Skip if chromophore2 is not one of chromophore1's neighbours
            # Also skip if chromophore2's ID is < chromophore1's ID to prevent
            # duplicates
            if (chromophore2.ID not in neighbours_ID) or (
                chromophore2.ID < chromophore1.ID
            ):
                continue
            # Update the orca input name
            input_name = chromophore1.orca_input.replace(
                ".inp", "-{:05d}.inp".format(chromophore2.ID)
            ).replace("single", "pair")
            # Find the correct relative image for the neighbour chromophore
            chromophore2_relative_image = neighbours_image[
                neighbours_ID.index(chromophore2.ID)
            ]
            chromophore2_transformation = list(
                np.array(chromophore1.image)
                - np.array(chromophore2.image)
                + np.array(chromophore2_relative_image)
            )
            # Find the dimer AAIDs and relative images for each atom
            AAIDs = chromophore1.AAIDs + chromophore2.AAIDs
            images = [[0, 0, 0] for i in range(len(chromophore1.AAIDs))] + [
                chromophore2_transformation for i in range(len(chromophore2.AAIDs))
            ]
            # Now add the terminating groups to both chromophores
            # Note that we would only ever expect both chromophores to require
            # termination or neither
            if chromophore1.terminate is True:
                term_group_posns1 = terminate_monomers(
                    chromophore1, parameter_dict, AA_morphology_dict
                )
                term_group_posns2 = terminate_monomers(
                    chromophore2, parameter_dict, AA_morphology_dict
                )
                # We don't want to add the terminating hydrogens for adjacent
                # monomers, so remove the ones that are within a particular
                # distance
                term_group_posns1, term_group_posns2 = remove_adjacent_terminators(
                    term_group_posns1, term_group_posns2
                )
                terminating_group_images1 = [
                    [0, 0, 0] for i in range(len(term_group_posns1))
                ]
                terminating_group_images2 = [
                    chromophore2_transformation for i in range(len(term_group_posns2))
                ]
                # Write the dimer input file
                write_orca_inp(
                    AA_morphology_dict,
                    AAIDs,
                    images,
                    term_group_posns1 + term_group_posns2,
                    terminating_group_images1 + terminating_group_images2,
                    "".join([parameter_dict["output_orca_directory"], input_name]),
                )
            else:
                # Write the dimer input file
                write_orca_inp(
                    AA_morphology_dict,
                    AAIDs,
                    images,
                    None,
                    None,
                    "".join([parameter_dict["output_orca_directory"], input_name]),
                )
    print("")


def remove_adjacent_terminators(group1, group2):
    pop_list = [[], []]
    for index1, terminating_hydrogen1 in enumerate(group1):
        for index2, terminating_hydrogen2 in enumerate(group2):
            separation = np.linalg.norm(terminating_hydrogen2 - terminating_hydrogen1)
            if separation < 1.2:
                pop_list[0].append(index1)
                pop_list[1].append(index2)
    for group_no, group in enumerate(pop_list):
        for index in sorted(group, reverse=True):
            try:
                [group1, group2][group_no].pop(index)
            except IndexError:
                raise SystemError(
                    "Tried to pop a termination group that does"
                    " not exist...are you sure this is an"
                    " atomistic morphology?"
                )
    return group1, group2


def write_orca_inp(
    AA_morphology_dict,
    AAIDs,
    images,
    terminating_group_posns,
    terminating_group_images,
    input_name,
):
    lines_to_write = []
    all_atom_types = []
    all_positions = []
    # Format the atom positions ready for orca
    for index, atom_ID in enumerate(AAIDs):
        # Cut the integer bit off the atomType. To allow atom types like Ca and
        # Br, where the atom is defined by one upper- and one lower-case letter,
        # use iter to return the first element of the list of upper case letters'
        # and the first element of the list of lower case letters in the atom
        # type (if any) and .join them together.
        atom_type = AA_morphology_dict["type"][atom_ID]
        all_atom_types.append(
            next(iter([char for char in atom_type if char.isupper()]), "")
            + next(iter([char for char in atom_type if char.islower()]), "")
        )
        # Add in the correct periodic images to the position
        all_positions.append(
            AA_morphology_dict["unwrapped_position"][atom_ID]
            + np.array(
                [
                    (
                        images[index][i]
                        * [
                            AA_morphology_dict["lx"],
                            AA_morphology_dict["ly"],
                            AA_morphology_dict["lz"],
                        ][i]
                    )
                    for i in range(3)
                ]
            )
        )
    # Now add in the terminating Hydrogens if necessary
    if terminating_group_posns is not None:
        for index, position in enumerate(terminating_group_posns):
            # Cut the integer bit off the atomType
            all_atom_types.append("H")
            # Add in the correct periodic images to the position
            all_positions.append(
                position
                + np.array(
                    [
                        (
                            terminating_group_images[index][i]
                            * [
                                AA_morphology_dict["lx"],
                                AA_morphology_dict["ly"],
                                AA_morphology_dict["lz"],
                            ][i]
                        )
                        for i in range(3)
                    ]
                )
            )
    # Now geometrically centralize all of the atoms that are to be included in
    # this input file to make it easier on orca
    central_position = np.array(
        [
            np.average(np.array(all_positions)[:, 0]),
            np.average(np.array(all_positions)[:, 1]),
            np.average(np.array(all_positions)[:, 2]),
        ]
    )
    # Create the lines to be written in the input file
    for index, position in enumerate(all_positions):
        lines_to_write.append(
            " {0:s}  {1:.5f}  {2:.5f}  {3:.5f}\n".format(
                all_atom_types[index],
                position[0] - central_position[0],
                position[1] - central_position[1],
                position[2] - central_position[2],
            )
        )
    # Load the orca input template
    orca_temp_dir = os.path.join(PROJECT_ROOT, "templates")
    orca_temp_test_dir = os.path.join(PROJECT_ROOT, "code", "unit_testing", "assets")
    try:
        with open(os.path.join(orca_temp_dir, "template.inp"), "r") as template_file:
            inp_file_lines = template_file.readlines()
    # In case running testbed:
    except FileNotFoundError:
        with open(
            os.path.join(orca_temp_test_dir, "template.inp"), "r"
        ) as template_file:
            inp_file_lines = template_file.readlines()
    # Insert the linesToWrite
    inp_file_lines[-1:-1] = lines_to_write
    # Write the orca input file
    with open(input_name, "w+") as orca_file:
        orca_file.writelines(inp_file_lines)
    print("\rOrca Input File written as", os.path.split(input_name[1]), end=" ")


def terminate_monomers(chromophore, parameter_dict, AA_morphology_dict):
    # No CG morphology, so we will use the UA -> AA code definition of which
    # atoms need to have hydrogens added to them.
    new_hydrogen_positions = []
    for atom_index_chromo, atom_index_morph in enumerate(chromophore.AAIDs):
        atom_type = AA_morphology_dict["type"][atom_index_morph]
        if atom_type not in parameter_dict["molecule_terminating_connections"].keys():
            continue
        bonded_AAIDs = []
        # Iterate over all termination connections defined for this atomType (in
        # case we are trying to do something mega complicated)
        for connection_info in parameter_dict["molecule_terminating_connections"][
            atom_type
        ]:
            for [bond_name, AAID1, AAID2] in chromophore.bonds:
                if AAID1 == atom_index_morph:
                    if AAID2 not in bonded_AAIDs:
                        bonded_AAIDs.append(AAID2)
                elif AAID2 == atom_index_morph:
                    if AAID1 not in bonded_AAIDs:
                        bonded_AAIDs.append(AAID1)
            if len(bonded_AAIDs) != connection_info[0]:
                continue
            new_hydrogen_positions += hf.get_terminating_positions(
                AA_morphology_dict["unwrapped_position"][atom_index_morph],
                [
                    AA_morphology_dict["unwrapped_position"][bonded_AAID]
                    for bonded_AAID in bonded_AAIDs
                ],
                1,
            )
    # Return terminatingGroups (positions of those hydrogens to be added to the
    # orca input)
    return new_hydrogen_positions


def get_orca_jobs(input_dir, parameter_dict, proc_IDs):
    # First delete any previous log files as we're about to start again with the
    # ZINDO/S calculations
    try:
        os.unlink(input_dir.replace("/input_orca", "/*.log"))
    except OSError:
        pass
    # Obtain a list of files to run
    single_orca_file_list = os.listdir(os.path.join(input_dir, "single"))
    pair_orca_file_list = os.listdir(os.path.join(input_dir, "pair"))
    orca_files_to_run = []
    for file_name in single_orca_file_list:
        if file_name[-4:] == ".inp":
            orca_files_to_run.append(os.path.join(input_dir, "single", file_name))
    for file_name in pair_orca_file_list:
        if file_name[-4:] == ".inp":
            orca_files_to_run.append(os.path.join(input_dir, "pair", file_name))
    orca_files_to_run.sort()
    if parameter_dict["overwrite_current_data"] is False:
        # Do not run any jobs that have already have an output file (and so have
        # at least started to
        # run if not finished)
        pop_list = []
        for job_no, job in enumerate(orca_files_to_run):
            try:
                with open(
                    job.replace("input_orca", "output_orca").replace(".inp", ".out"),
                    "r",
                ):
                    pop_list.append(job_no)
            except IOError:
                pass
        pop_list.sort(reverse=True)
        for pop_index in pop_list:
            orca_files_to_run.pop(pop_index)
    if len(orca_files_to_run) == 0:
        return []
    # Create a jobslist for each procID
    jobs_list = [
        orca_files_to_run[
            i : i + (int(np.ceil(len(orca_files_to_run) / len(proc_IDs))))
        ]
        for i in range(
            0,
            len(orca_files_to_run),
            int(np.ceil(len(orca_files_to_run) / float(len(proc_IDs)))),
        )
    ]
    return jobs_list


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
    create_input_files(chromophore_list, AA_morphology_dict, parameter_dict)
    input_dir = os.path.join(
        parameter_dict["output_orca_directory"], "chromophores", "input_orca"
    )
    proc_IDs = parameter_dict["proc_IDs"]
    jobs_list = get_orca_jobs(input_dir, parameter_dict, proc_IDs)
    # Shuffle the jobsList to spread it out over the cores
    np.random.shuffle(jobs_list)
    number_of_inputs = sum([len(orca_files_to_run) for orca_files_to_run in jobs_list])
    print("Found", number_of_inputs, "orca files to run.")
    if number_of_inputs > 0:
        # Create pickle file containing the jobs sorted by ProcID to be picked
        # up by single_core_run_orca.py
        pickle_name = input_dir.replace("input_orca", "orca_jobs.pickle")
        with open(pickle_name, "wb+") as pickle_file:
            pickle.dump(jobs_list, pickle_file)
        print("Orca jobs list written to", pickle_name)
        if len(jobs_list) <= len(proc_IDs):
            proc_IDs = proc_IDs[: len(jobs_list)]
        running_jobs = []
        # Open the required processes to execute the orca jobs
        for CPU_rank, jobs in enumerate(jobs_list):
            running_jobs.append(
                sp.Popen(
                    [
                        "python",
                        SINGLE_ORCA_RUN_FILE,
                        parameter_dict["output_orca_directory"],
                        parameter_dict["output_morphology_directory"],
                        str(CPU_rank),
                        str(int(parameter_dict["overwrite_current_data"])),
                        str(int(parameter_dict["remove_orca_inputs"])),
                    ]
                )
            )
        # Wait for all jobs to complete
        [p.wait() for p in running_jobs]
        # Delete the job pickle
        os.system(" ".join(["rm", pickle_name]))
    return (
        AA_morphology_dict,
        CG_morphology_dict,
        CG_to_AAID_master,
        parameter_dict,
        chromophore_list,
    )


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
