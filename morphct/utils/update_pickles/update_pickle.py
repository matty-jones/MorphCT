import re
from morphct.code import helper_functions as hf
from morphct.code import obtain_chromophores as oc
from morphct.templates import par_template


def convert_params(old_parameter_dict):
    # First, rename all of the variables to the new format
    renamed_parameter_dict = add_underscores(old_parameter_dict)
    # Then, reorganise any redundant parameters (subspecies stuff)
    fixed_parameter_dict = rename_old(renamed_parameter_dict)
    # Remove any redundant parameters from the current dict
    modified_parameter_dict = remove_redundant(fixed_parameter_dict)
    # Finally, add any missing parameters from the template
    added_parameter_dict = add_missing_parameters(modified_parameter_dict)
    # Rewrite the parameter file?
    # Return the parameter dictionary to be repickled
    exit()


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
    print("In par dict but not template =", response - expected)
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
        old_parameter_dict['input_morph_file'] = old_parameter_dict.pop('input_morphology_file')
    except KeyError:
        pass
    try:
        old_parameter_dict['output_morph_file'] = old_parameter_dict.pop('output_morphology_file')
    except KeyError:
        pass
    # Do subspecies
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
    print("In template but not par dict =", missing_params)
    for key in missing_params:
        value = eval('.'.join(['par_template', key]))
        old_parameter_dict[key] = value
    exit()


if __name__ == "__main__":
    old_pickle_file = './acceptorCrystal.pickle'
    pickle_data = hf.load_pickle(old_pickle_file)
    AA_morphology_dict = pickle_data[0]
    CG_morphology_dict = pickle_data[1]
    CG_to_AAID_master = pickle_data[2]
    old_parameter_dict = pickle_data[3]
    old_chromophore_list = pickle_data[4]

    new_parameter_dict = convert_params(old_parameter_dict)
