import argparse
import os
import pickle
import re
import shutil
import sys
import glob
import numpy as np
from morphct.code import helper_functions as hf
from morphct.code import obtain_chromophores as oc
from morphct.definitions import PROJECT_ROOT
from morphct.templates import par as par_template


def convert_params(old_parameter_dict, file_name):
    # First, rename all of the variables to the new format
    print("Updating parameter names to PEP8 format...")
    new_parameter_dict = add_underscores(old_parameter_dict)
    # Then, reorganise any redundant parameters (subspecies stuff)
    print("Reorganising parameters...")
    new_parameter_dict = rename_old(new_parameter_dict)
    # Remove any redundant parameters from the current dict
    print("Removing deprecated parameters...")
    new_parameter_dict = remove_redundant(new_parameter_dict)
    # Finally, add any missing parameters from the template
    print("Setting missing parameters to defaults...")
    new_parameter_dict = add_missing_parameters(new_parameter_dict)
    # Rewrite the parameter file
    if file_name is not None:
        print("Rewriting parameter file...")
        rewrite_parameter_file(new_parameter_dict, file_name)
    # Return the parameter dictionary to be repickled
    return new_parameter_dict


def add_underscores(old_parameter_dict):
    leave_these_capitalised = [
        "MA",
        "CG",
        "ID",
        "AA",
        "TI",
        "HOMO",
        "LUMO",
        "KMC",
        "VRH",
        "MO",
        "DOS",
        "AAID",
    ]
    leave_these_capitalised += [
        "".join([abbrev, "s"]) for abbrev in leave_these_capitalised
    ]
    new_parameter_dict = {}
    for key, value in old_parameter_dict.items():
        # Some hardcoded ones because I can't work out the regex
        # to detect them nicely
        if "DoSSTD" in key:
            new_key = "_DOS_std_".join(key.split("DoSSTD"))
        else:
            catch_lower = re.compile("([a-z])([A-Z]+)").sub(r"\1_\2", key)
            catch_upper = re.compile("([A-Z])([A-Z](?![s])[a-z])").sub(
                r"\1_\2", catch_lower
            )
            split_key = catch_upper.split("_")
            new_key = "_".join(
                [
                    el.lower() if el not in leave_these_capitalised else el
                    for el in split_key
                ]
            )
        new_parameter_dict[new_key] = value
    return new_parameter_dict


def rename_old(old_parameter_dict):
    expected = set(
        [parameter for parameter in dir(par_template) if parameter[0] != "_"]
    )
    response = set(old_parameter_dict.keys())
    # Make following changes:
    #   input_morphology_directory -> input_morph_dir
    #   output_morphology_directory -> output_morph_dir
    #   execute_zindo -> execute_ZINDO
    #   execute_finegraining -> execute_fine_graining
    #   input_morphology_file -> input_morph_file
    #   output_morphology_file -> output_morph_file
    try:
        old_parameter_dict["input_morph_dir"] = old_parameter_dict.pop(
            "input_morphology_directory"
        )
    except KeyError:
        pass
    try:
        old_parameter_dict["output_morph_dir"] = old_parameter_dict.pop(
            "output_morphology_directory"
        )
    except KeyError:
        pass
    try:
        old_parameter_dict["execute_fine_graining"] = old_parameter_dict.pop(
            "execute_finegraining"
        )
    except KeyError:
        pass
    try:
        old_parameter_dict["execute_ZINDO"] = old_parameter_dict.pop("execute_zindo")
    except KeyError:
        pass
    try:
        old_parameter_dict["morphology"] = os.path.split(
            old_parameter_dict.pop("input_morphology_file")
        )[1]
    except KeyError:
        pass
    try:
        old_parameter_dict["output_morph_file"] = old_parameter_dict.pop(
            "output_morphology_file"
        )
    except KeyError:
        pass
    # Do subspecies
    # Check if this is from when we only had one reorganisation_energy
    try:
        old_parameter_dict["reorganisation_energy_donor"]
    except KeyError:
        old_parameter_dict["reorganisation_energy_donor"] = old_parameter_dict[
            "reorganisation_energy"
        ]
        old_parameter_dict["reorganisation_energy_acceptor"] = 0.0000
    # Check the CG_site_species to see if the chromophore species are
    # capitalised or not
    # if 'Donor' in old_parameter_dict['CG_site_species'].keys():
    #    donor_species = 'Donor'
    # else:
    #    donor_species = 'donor'
    # if 'Acceptor' in old_parameter_dict['CG_site_species'].keys():
    #    acceptor_species = 'Acceptor'
    # else:
    #    acceptor_species = 'acceptor'
    try:
        chromophore_species = {
            "Donor": {
                "target_DOS_std": old_parameter_dict.pop("target_DOS_std_HOMO"),
                "literature_MO": old_parameter_dict.pop("literature_HOMO"),
                "VRH_delocalisation": 2e-10,
                "species": "donor",  # donor_species,
                "reorganisation_energy": old_parameter_dict.pop(
                    "reorganisation_energy_donor"
                ),
            },
            "Acceptor": {
                "target_DOS_std": old_parameter_dict.pop("target_DOS_std_LUMO"),
                "literature_MO": old_parameter_dict.pop("literature_LUMO"),
                "VRH_delocalisation": 2e-10,
                "species": "acceptor",  # acceptor_species,
                "reorganisation_energy": old_parameter_dict.pop(
                    "reorganisation_energy_acceptor"
                ),
            },
        }
        old_parameter_dict["chromophore_species"] = chromophore_species
    except KeyError:
        pass
    expected = set(
        [parameter for parameter in dir(par_template) if parameter[0] != "_"]
    )
    response = set(old_parameter_dict.keys())
    return old_parameter_dict


def remove_redundant(old_parameter_dict):
    remove_these_keys = [
        "electrical_field",
        "execute_extract_molecules",
        "koopmans_hopping_prefactor",
    ]
    for key in remove_these_keys:
        try:
            old_parameter_dict.pop(key)
        except KeyError:
            pass
    return old_parameter_dict


def add_missing_parameters(old_parameter_dict):
    expected = set(
        [parameter for parameter in dir(par_template) if parameter[0] != "_"]
    )
    response = set(old_parameter_dict.keys())
    missing_params = expected - response
    for key in missing_params:
        value = eval(".".join(["par_template", key]))
        old_parameter_dict[key] = value
    return old_parameter_dict


def rewrite_parameter_file(new_parameter_dict, new_parameter_file):
    # Then update the empty template with all of the right variables
    with open(
        os.path.join(os.path.split(par_template.__file__)[0], "empty_par_template.py"),
        "r",
    ) as empty_file:
        lines = empty_file.readlines()
    for line_number, line in enumerate(lines):
        if " =\n" not in line:
            continue
        split_line = line.split(" =")
        split_line.insert(1, repr(new_parameter_dict[split_line[0]]))
        split_line.insert(1, " = ")
        lines[line_number] = "".join(split_line)
    # Now write that file
    with open(new_parameter_file, "w+") as par_file:
        par_file.writelines(lines)
    print("Updated parameter written to", new_parameter_file)


def convert_chromos(
    old_chromophore_list,
    CG_morphology_dict,
    AA_morphology_dict,
    CG_to_AAID_master,
    parameter_dict,
    sim_dims,
):
    new_chromophore_list = []
    print("Updating the chromophore list (this might take a few minutes)...")
    for old_chromo in old_chromophore_list:
        print(
            "\rUpdating chromophore {0:d} of {1:d}...".format(
                old_chromo.ID + 1, len(old_chromophore_list)
            ),
            end=" ",
        )
        # Set up an empty chromophore instance using the parameter_dict
        new_chromo = oc.chromophore(
            old_chromo.ID,
            old_chromo.CGIDs,
            CG_morphology_dict,
            AA_morphology_dict,
            CG_to_AAID_master,
            parameter_dict,
            sim_dims,
        )
        # Copy over the old_chromophore properties for the following:
        # super_cell_positions, super_cell_images
        # Energy levels (HOMO_1, HOMO, LUMO, LUMO_1)
        # Neighbours (neighbours, dissociation_neighbours, neighbours_TI,
        # neighbours_delta_E)
        new_chromo = update_new_chromo(old_chromo, new_chromo, sim_dims)
        new_chromophore_list.append(new_chromo)
    print()
    return new_chromophore_list


def update_new_chromo(old_chromo, new_chromo, sim_dims):
    # Create a properties dict that maps from the old parameter to the new one
    properties = {
        "superCellPosns": "super_cell_posns",
        "superCellImages": "super_cell_images",
        "HOMO_1": "HOMO_1",
        "HOMO": "HOMO",
        "LUMO": "LUMO",
        "LUMO_1": "LUMO_1",
        "neighbours": "neighbours",
        "dissociationNeighbours": "dissociation_neighbours",
        "neighboursTI": "neighbours_TI",
        "neighboursDeltaE": "neighbours_delta_E",
    }
    if ("superCellPosns" not in old_chromo.__dict__) or (
        "superCellImages" not in old_chromo.__dict__
    ):
        old_chromo = add_super_cell_data(old_chromo, sim_dims)
    if "dissociationNeighbours" not in old_chromo.__dict__:
        old_chromo.__dict__["dissociationNeighbours"] = []
    for old_prop, new_prop in properties.items():
        new_chromo.__dict__[new_prop] = old_chromo.__dict__[old_prop]
    return new_chromo


def add_super_cell_data(chromophore, sim_dims):
    box = np.array([axis[1] - axis[0] for axis in sim_dims])
    chromophore.superCellImages = [
        np.array([x, y, z])
        for x in range(-1, 2)
        for y in range(-1, 2)
        for z in range(-1, 2)
    ]
    chromophore.superCellPosns = [
        np.array(chromophore.posn) + (box * image)
        for image in chromophore.superCellImages
    ]
    return chromophore


def get_parameter_file(directory):
    parameter_files = []
    directory_files = [
        file_name
        for file_name in os.listdir(directory)
        if os.path.isfile(os.path.join(directory, file_name))
    ]
    for file_name in directory_files:
        try:
            with open(os.path.join(directory, file_name), "r") as file_handle:
                for line in file_handle:
                    if "MorphCT.simulation" in line:
                        parameter_files.append(file_name)
        except UnicodeDecodeError:
            continue
    # Remove previous backups
    parameter_files = [
        file_name for file_name in parameter_files if ".bak_par" not in file_name
    ]
    if len(parameter_files) == 0:
        print("No parameter file found to back up.")
        return None
    elif len(parameter_files) > 1:
        print(
            "Found",
            len(parameter_files),
            "appropriate parameter files:",
            parameter_files,
        )
        print("Backing up all of them just in case...")
    return parameter_files


def load_pickle_data(old_pickle_file, directory):
    try:
        pickle_data = hf.load_pickle(old_pickle_file)
    except SyntaxError:
        # This pickle is from a python 2 version of MorphCT. We can inject
        # a line into the obtainChromophores.py to try and fix the import
        # error
        print("Pickle is from Python 2, updating the code so it can be imported...")
        code_file_names = ["obtainChromophores.py", "helperFunctions.py"]
        for code_file_name in code_file_names:
            print("Updating", "".join([code_file_name, "..."]))
            with open(os.path.join(directory, code_file_name), "r") as code_file:
                code_lines = code_file.readlines()
            for line_number, line in enumerate(code_lines):
                if "import cPickle as pickle" in line:
                    code_lines[line_number] = "import pickle\n"
                elif ("#" not in line) and ("print" in line):
                    n_spaces = len(line) - len(line.lstrip())
                    code_lines[line_number] = "".join([" " * n_spaces, "print()\n"])
            with open(os.path.join(directory, code_file_name), "w+") as code_file:
                code_file.writelines(code_lines)
        pickle_data = hf.load_pickle(old_pickle_file)
    return pickle_data


def convert_KMC(morphology_directory):
    KMC_data = load_KMC_results_pickle(morphology_directory)
    print("Updating KMC result keys to PEP8 format...")
    KMC_data = add_underscores(KMC_data)
    print("Rewriting KMC results file...")
    new_file_name = os.path.join(morphology_directory, "KMC", "KMC_results.pickle")
    with open(new_file_name, "wb+") as pickle_file:
        pickle.dump(KMC_data, pickle_file)
    print("Updated KMC results file written to", new_file_name)


def load_KMC_results_pickle(directory):
    KMC_pickle = os.path.join(directory, "KMC", "KMCResults.pickle")
    try:
        with open(KMC_pickle, "rb") as pickle_file:
            carrier_data = pickle.load(pickle_file)
    except FileNotFoundError:
        print("No final KMC_results.pickle found. Creating it from incomplete parts...")
        create_results_pickle(directory)
        with open(KMC_pickle, "rb") as pickle_file:
            carrier_data = pickle.load(pickle_file)
    except UnicodeDecodeError:
        with open(KMC_pickle, "rb") as pickle_file:
            carrier_data = pickle.load(pickle_file, encoding="latin1")
    return carrier_data


def create_results_pickle(directory):
    cores_list = []
    for file_name in glob.glob(os.path.join(directory, "KMC", "*")):
        try:
            cores_list.append(re.findall("([_])(..)([\.])", file_name)[0][1])
        except IndexError:
            pass
    cores_list = sorted(list(set(cores_list)))
    results_pickles_list = []
    keep_list = []
    for core in cores_list:
        # Check if there is already a finished KMC_results pickle
        main = os.path.join(
            directory, "KMC", "KMCResults_{:02d}.pickle".format(int(core))
        )
        if os.path.exists(main):
            results_pickles_list.append(main)
            keep_list.append(None)
            continue
        # If not, find the slot1 and slot2 pickle that is most recent
        slot1 = os.path.join(
            directory, "KMC", "KMCSlot1Results_{:02d}.pickle".format(int(core))
        )
        slot2 = os.path.join(
            directory, "KMC", "KMCSlot2Results_{:02d}.pickle".format(int(core))
        )
        if os.path.exists(slot1) and not os.path.exists(slot2):
            keep_list.append(slot1)
        elif os.path.exists(slot2) and not os.path.exists(slot1):
            keep_list.append(slot2)
        elif os.path.getsize(slot1) >= os.path.getsize(slot2):
            keep_list.append(slot1)
        else:
            keep_list.append(slot2)
    print("{:d} pickle files found to combine!".format(len(keep_list)))
    print("Combining", keep_list)
    for keeper in zip(cores_list, keep_list):
        # Skip this core if we already have a finished KMC_results for it
        if keeper[1] is None:
            continue
        new_name = os.path.join(
            directory, "KMC", "KMC_results_{}.pickle".format(keeper[0])
        )
        shutil.copyfile(str(keeper[1]), new_name)
        results_pickles_list.append(new_name)
    combine_results_pickles(directory, results_pickles_list)


def combine_results_pickles(directory, pickle_files):
    combined_data = {}
    pickle_files = sorted(pickle_files)
    for file_name in pickle_files:
        # The pickle was repeatedly dumped to, in order to save time.
        # Each dump stream is self-contained, so iteratively unpickle to add the new data.
        with open(file_name, "rb") as pickle_file:
            pickled_data = pickle.load(pickle_file)
            for key, val in pickled_data.items():
                if val is None:
                    continue
                if key not in combined_data:
                    combined_data[key] = val
                else:
                    combined_data[key] += val
    # Write out the combined data
    print("Writing out the combined pickle file...")
    combined_file_loc = os.path.join(directory, "KMC", "KMC_results.pickle")
    with open(combined_file_loc, "wb+") as pickle_file:
        pickle.dump(combined_data, pickle_file)
    print("Complete data written to", combined_file_loc)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p",
        "--pickle_mode",
        action="store_true",
        required=False,
        help="""When this flag is passed, `pickle mode' is activated
                        in update_pickles. In this mode, the user does not supply a
                        morphology or morphology/code directory to operate on, but
                        instead provides a pickle file.
                        Note that the correct version of obtainChromophores.py and
                        helperFunctions.py must also be present within the same
                        directory as the pickle.
                        In this mode, the pickle file is still backed up, but the
                        parameter file is not written out.""",
    )
    args, input_list = parser.parse_known_args()
    if args.pickle_mode:
        print("Update_pickle is running in `pickle mode' (see -h for more details).")
    else:
        print("Update_pickle is running in `directory mode' (default).")
    # Iterate over all input strings (directories or pickles)
    for input_string in input_list:
        print()
        if args.pickle_mode:
            # Load the pickle directly, and do not output a parameter file
            # NOTE: It looks like old_ and new_ are the wrong way around here
            # and below, but it will be correct after the shutil.copy has taken
            # place (we want to move the old pickle to .bak_pickle, and then
            # overwrite the .pickle with the new data)
            old_pickle_file = input_string.replace(".pickle", ".bak_pickle")
            new_pickle_file = input_string
            directory = os.path.dirname(os.path.abspath(input_string))
            new_parameter_file = None
        else:
            # Make sure that we're in the code directory of the input morphology
            if len(input_string.split("/code")) == 1:
                directory = os.path.join(input_string, "code")
            else:
                directory = input_string
            morphology_name = os.path.dirname(directory).split("/")[-1]
            old_pickle_file = os.path.join(
                directory, "".join([morphology_name, ".bak_pickle"])
            )
            new_pickle_file = os.path.join(
                directory, "".join([morphology_name, ".pickle"])
            )
            # Find the parameter file and back it up
            parameter_files = get_parameter_file(directory)
            if (parameter_files is not None) and (len(parameter_files) > 0):
                parameter_found = False
                for parameter_file in parameter_files:
                    # Make a copy of the old parameter_dict
                    if parameter_found is False:
                        old_parameter_file = os.path.join(
                            directory, parameter_file.replace(".py", ".bak_par")
                        )
                        new_parameter_file = os.path.join(directory, parameter_file)
                        print("Found parameter file at", new_parameter_file)
                        parameter_found = True
                    if os.path.isfile(old_parameter_file):
                        print(
                            "Backup parameter file already exists at",
                            old_parameter_file,
                        )
                        print("Skipping creating a new backup...")
                    else:
                        print(
                            "Backing up", new_parameter_file, "to", old_parameter_file
                        )
                        shutil.copy(new_parameter_file, old_parameter_file)
            else:
                new_parameter_file = os.path.join(
                    directory, "".join(["par_", morphology_name, ".py"])
                )
        # Now back up the pickle before we overwrite it with the new data
        # (both pickle_mode and dir_mode)
        print("Considering morphology", directory)
        try:
            print("Found pickle at", new_pickle_file)
            # Make a copy of the old pickle
            if os.path.isfile(old_pickle_file):
                print("Backup pickle already exists at", old_pickle_file)
                print("Skipping creating a new backup...")
            else:
                print("Backing up", new_pickle_file, "to", old_pickle_file)
                shutil.copy(new_pickle_file, old_pickle_file)
        except FileNotFoundError:
            # This exception can never be raised in pickle_mode (so we don't
            # need to worry about morphology_name not being defined)
            error_string = "".join(
                [
                    "Tried to find ",
                    morphology_name,
                    ".pickle in ",
                    directory,
                    " but couldn't.",
                ]
            )
            raise SystemError(error_string)
        # Use the current MorphCT stored in the current directory
        sys.path.insert(0, os.path.abspath(directory))
        try:
            pickle_data = load_pickle_data(old_pickle_file, directory)
        except ImportError:
            # Trying to load a pickle from the MorphCT package pre 3.0
            print(
                "Importing pickle has failed, most likely a pickle from package-style MorphCT from before"
                " the PEP8 changes (v3.0.0). Copying code to help unpickling..."
            )
            shutil.copy(
                os.path.join(directory, "helperFunctions.py"),
                os.path.join(PROJECT_ROOT, "code", "helperFunctions.py"),
            )
            shutil.copy(
                os.path.join(directory, "obtainChromophores.py"),
                os.path.join(PROJECT_ROOT, "code", "obtainChromophores.py"),
            )
            pickle_data = load_pickle_data(old_pickle_file, directory)
            print("Deleting copied code to keep things tidy...")
            os.remove(os.path.join(PROJECT_ROOT, "code", "helperFunctions.py"))
            os.remove(os.path.join(PROJECT_ROOT, "code", "obtainChromophores.py"))
        AA_morphology_dict = pickle_data[0]
        CG_morphology_dict = pickle_data[1]
        CG_to_AAID_master = pickle_data[2]
        old_parameter_dict = pickle_data[3]
        old_chromophore_list = pickle_data[4]
        sim_dims = [
            [-AA_morphology_dict["lx"] / 2.0, AA_morphology_dict["lx"] / 2.0],
            [-AA_morphology_dict["ly"] / 2.0, AA_morphology_dict["ly"] / 2.0],
            [-AA_morphology_dict["lz"] / 2.0, AA_morphology_dict["lz"] / 2.0],
        ]

        # Update the parameter dict and pickle files to include the new data
        new_parameter_dict = convert_params(old_parameter_dict, new_parameter_file)

        # #DEBUG
        # # Temp for Mike jobs
        # new_parameter_dict['AA_rigid_body_species'] = {'Acceptor': list(range(0, 1970, 2)),
        #                                                'Donor': list(range(1, 1970, 2))}
        # print("DANGER: INJECTED HARDCODE TO MAKE BDT-TPD WORK")
        # #

        new_chromophore_list = convert_chromos(
            old_chromophore_list,
            CG_morphology_dict,
            AA_morphology_dict,
            CG_to_AAID_master,
            new_parameter_dict,
            sim_dims,
        )
        # Write out the new data
        hf.write_pickle(
            [
                AA_morphology_dict,
                CG_morphology_dict,
                CG_to_AAID_master,
                new_parameter_dict,
                new_chromophore_list,
            ],
            new_pickle_file,
        )

        # Update the KMC data
        if args.pickle_mode is False:
            morphology_directory = directory.replace("/code", "")
            convert_KMC(morphology_directory)

        # Remove the current directory from the path, ready for the next
        # directory
        sys.path.pop(0)


if __name__ == "__main__":
    main()
