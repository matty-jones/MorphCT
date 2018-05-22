import argparse
import os
import re
import shutil
import sys
import numpy as np
from morphct.code import helper_functions as hf
from morphct.code import obtain_chromophores as oc
from morphct.templates import par_template


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
    print(new_parameter_dict.keys())
    exit()
    print("Rewriting parameter file...")
    rewrite_parameter_file(new_parameter_dict, file_name)
    # Return the parameter dictionary to be repickled
    return new_parameter_dict


def add_underscores(old_parameter_dict):
    leave_these_capitalised = ['MA', 'CG', 'ID', 'AA', 'TI', 'HOMO', 'LUMO', 'KMC', 'VRH', 'MO', 'DOS', 'AAID']
    leave_these_capitalised += [abbrev + 's' for abbrev in leave_these_capitalised]
    new_parameter_dict = {}
    for key, value in old_parameter_dict.items():
        # Some hardcoded ones because I can't work out the regex
        # to detect them nicely
        if 'DoSSTD' in key:
            new_key = '_DOS_std_'.join(key.split('DoSSTD'))
        else:
            catch_lower = re.compile("([a-z])([A-Z]+)").sub(r"\1_\2", key)
            catch_upper = re.compile("([A-Z])([A-Z](?![s])[a-z])").sub(r"\1_\2", catch_lower)
            split_key = catch_upper.split('_')
            new_key = '_'.join([el.lower() if el not in leave_these_capitalised else el for el in split_key])
        new_parameter_dict[new_key] = value
    return new_parameter_dict


def rename_old(old_parameter_dict):
    expected = set([parameter for parameter in dir(par_template) if parameter[0] != '_'])
    response = set(old_parameter_dict.keys())
    # Make following changes:
    #   input_morphology_directory -> input_morph_dir
    #   output_morphology_directory -> output_morph_dir
    #   execute_zindo -> execute_ZINDO
    #   execute_finegraining -> execute_fine_graining
    #   input_morphology_file -> input_morph_file
    #   output_morphology_file -> output_morph_file
    try:
        old_parameter_dict['input_morph_dir'] = old_parameter_dict.pop('input_morphology_directory')
    except KeyError:
        pass
    try:
        old_parameter_dict['output_morph_dir'] = old_parameter_dict.pop('output_morphology_directory')
    except KeyError:
        pass
    try:
        old_parameter_dict['execute_fine_graining'] = old_parameter_dict.pop('execute_finegraining')
    except KeyError:
        pass
    try:
        old_parameter_dict['execute_ZINDO'] = old_parameter_dict.pop('execute_zindo')
    except KeyError:
        pass
    try:
        old_parameter_dict['morphology'] = os.path.split(old_parameter_dict.pop('input_morphology_file'))[1]
    except KeyError:
        pass
    try:
        old_parameter_dict['output_morph_file'] = old_parameter_dict.pop('output_morphology_file')
    except KeyError:
        pass
    # Do subspecies
    try:
        chromophore_species = {'Donor':
                               {'target_DOS_std': old_parameter_dict.pop('target_DOS_std_HOMO'),
                                'literature_MO': old_parameter_dict.pop('literature_HOMO'),
                                'VRH_delocalisation': 2e-10,
                                'species': 'donor',
                                'reorganisation_energy': old_parameter_dict.pop('reorganisation_energy_donor')
                               },
                               'Acceptor':
                               {'target_DOS_std': old_parameter_dict.pop('target_DOS_std_LUMO'),
                                'literature_MO': old_parameter_dict.pop('literature_LUMO'),
                                'VRH_delocalisation': 2e-10,
                                'species': 'acceptor',
                                'reorganisation_energy': old_parameter_dict.pop('reorganisation_energy_acceptor')
                               }
                              }
        old_parameter_dict['chromophore_species'] = chromophore_species
    except KeyError:
        pass
    expected = set([parameter for parameter in dir(par_template) if parameter[0] != '_'])
    response = set(old_parameter_dict.keys())
    return old_parameter_dict


def remove_redundant(old_parameter_dict):
    remove_these_keys = ['electrical_field', 'execute_extract_molecules']
    for key in remove_these_keys:
        try:
            old_parameter_dict.pop(key)
        except KeyError:
            pass
    return old_parameter_dict


def add_missing_parameters(old_parameter_dict):
    expected = set([parameter for parameter in dir(par_template) if parameter[0] != '_'])
    response = set(old_parameter_dict.keys())
    missing_params = expected - response
    for key in missing_params:
        value = eval('.'.join(['par_template', key]))
        old_parameter_dict[key] = value
    return old_parameter_dict


def rewrite_parameter_file(new_parameter_dict, new_parameter_file):
    # Then update the empty template with all of the right variables
    with open(os.path.join(os.path.split(par_template.__file__)[0], 'empty_par_template.py'), 'r') as empty_file:
        lines = empty_file.readlines()
    for line_number, line in enumerate(lines):
        if ' =\n' not in line:
            continue
        split_line = line.split(' =')
        split_line.insert(1, repr(new_parameter_dict[split_line[0]]))
        split_line.insert(1, ' = ')
        lines[line_number] = ''.join(split_line)
    # Now write that file
    with open(new_parameter_file, 'w+') as par_file:
        par_file.writelines(lines)
    print("Updated parameter written to", new_parameter_file)


def convert_chromos(old_chromophore_list, CG_morphology_dict, AA_morphology_dict, CG_to_AAID_master, parameter_dict,
                    sim_dims):
    new_chromophore_list = []
    print("Updating the chromophore list (this might take a few minutes)...")
    for old_chromo in old_chromophore_list:
        print("\rUpdating chromophore", old_chromo.ID, end=' ')
        # Set up an empty chromophore instance using the parameter_dict
        new_chromo = oc.chromophore(old_chromo.ID, old_chromo.CGIDs, CG_morphology_dict, AA_morphology_dict,
                                    CG_to_AAID_master, parameter_dict, sim_dims)
        # Copy over the old_chromophore properties for the following:
        # super_cell_positions, super_cell_images
        # Energy levels (HOMO_1, HOMO, LUMO, LUMO_1)
        # Neighbours (neighbours, dissociation_neighbours, neighbours_TI,
        # neighbours_delta_E)
        new_chromo = update_new_chromo(old_chromo, new_chromo)
        new_chromophore_list.append(new_chromo)
    return new_chromophore_list


def update_new_chromo(old_chromo, new_chromo):
    # Create a properties dict that maps from the old parameter to the new one
    properties = {'superCellPosns': 'super_cell_posns',
                  'superCellImages': 'super_cell_images',
                  'HOMO_1': 'HOMO_1',
                  'HOMO': 'HOMO',
                  'LUMO': 'LUMO',
                  'LUMO_1': 'LUMO_1',
                  'neighbours': 'neighbours',
                  'dissociationNeighbours': 'dissociation_neighbours',
                  'neighboursTI': 'neighbours_TI',
                  'neighboursDeltaE': 'neighbours_delta_E',
                 }
    for old_prop, new_prop in properties.items():
        new_chromo.__dict__[new_prop] = old_chromo.__dict__[old_prop]
    return new_chromo


def get_parameter_file(directory):
    parameter_files = []
    directory_files = [file_name for file_name in os.listdir(directory) if
                       os.path.isfile(os.path.join(directory, file_name))]
    for file_name in directory_files:
        try:
            with open(os.path.join(directory, file_name), 'r') as file_handle:
                for line in file_handle:
                    if 'MorphCT.simulation' in line:
                        parameter_files.append(file_name)
        except UnicodeDecodeError:
            continue
    # Remove previous backups
    parameter_files = [file_name for file_name in parameter_files if '.bak_par' not in file_name]
    if len(parameter_files) == 0:
        print("No parameter file found to back up.")
        return None
    elif len(parameter_files) > 1:
        print("Found", len(parameter_files), "appropriate parameter files:", parameter_files)
        print("Backing up all of them just in case...")
    return parameter_files


def main():
    parser = argparse.ArgumentParser()
    args, directory_list = parser.parse_known_args()
    for raw_dir in directory_list:
        # Make sure that we're in the code directory of the input morphology
        if len(raw_dir.split('/code')) == 1:
            directory = os.path.join(raw_dir, 'code')
        else:
            directory = raw_dir
        # Use the current MorphCT stored in the current directory
        sys.path.insert(0, os.path.abspath(directory))
        morphology_name = os.path.dirname(directory).split('/')[-1]
        print("Considering morphology directory", directory)
        try:
            # NOTE: It looks like old_ and new_ are the wrong way around here,
            # but it will be correct after the shutil.copy has taken place
            # (we want to move the old pickle to .bak_pickle, and then
            # overwrite the .pickle with the new data)
            old_pickle_file = os.path.join(directory, ''.join([morphology_name, '.bak_pickle']))
            new_pickle_file = os.path.join(directory, ''.join([morphology_name, '.pickle']))
            print("Found pickle at", new_pickle_file)
            # Make a copy of the old pickle
            print("Backing up", new_pickle_file, "to", old_pickle_file)
            shutil.copy(new_pickle_file, old_pickle_file)
        except:
            error_string = ''.join(['Tried to find ', morphology_name, '.pickle in ', directory, " but couldn't."])
            raise SystemError(error_string)
        parameter_files = get_parameter_file(directory)
        if parameter_files is not None:
            parameter_found = False
            for parameter_file in parameter_files:
                # Make a copy of the old parameter_dict
                if parameter_found is False:
                    old_parameter_file = os.path.join(directory, parameter_file.replace('.py', '.bak_par'))
                    new_parameter_file = os.path.join(directory, parameter_file)
                    print("Found parameter file at", new_parameter_file)
                    parameter_found = True
                shutil.copy(os.path.join(directory, parameter_file), os.path.join(
                    directory, parameter_file.replace('.py', '.bak_par')))
        pickle_data = hf.load_pickle(old_pickle_file)
        AA_morphology_dict = pickle_data[0]
        CG_morphology_dict = pickle_data[1]
        CG_to_AAID_master = pickle_data[2]
        old_parameter_dict = pickle_data[3]
        old_chromophore_list = pickle_data[4]
        sim_dims = [[-AA_morphology_dict['lx'] / 2.0, AA_morphology_dict['lx'] / 2.0],
                    [-AA_morphology_dict['ly'] / 2.0, AA_morphology_dict['ly'] / 2.0],
                    [-AA_morphology_dict['lz'] / 2.0, AA_morphology_dict['lz'] / 2.0]]
        # Update the parameter dict and pickle files to include the new data
        new_parameter_dict = convert_params(old_parameter_dict, new_parameter_file)
        new_chromophore_list = convert_chromos(old_chromophore_list, CG_morphology_dict, AA_morphology_dict,
                                               CG_to_AAID_master, new_parameter_dict, sim_dims)
        # Write out the new data
        hf.write_pickle([AA_morphology_dict, CG_morphology_dict, CG_to_AAID_master, new_parameter_dict,
                         new_chromophore_list],
                        new_pickle_file)
        # Remove the current directory from the path, ready for the next
        # directory
        sys.path.pop(0)


if __name__ == "__main__":
    main()
