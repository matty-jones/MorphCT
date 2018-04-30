from morphct.code import obtain_chromophores
from morphct.code import execute_ZINDO
from morphct.code import helper_functions as hf
import copy
import random as R
import numpy as np


def test_find_neighbours(pickle_file):
    pickle_data = hf.load_pickle(pickle_file)
    AA_morphology_dict = pickle_data[0]
    parameter_dict = pickle_data[3]
    old_chromophore_list = pickle_data[4]
    empty_cut_off_chromophore_list = copy.deepcopy(old_chromophore_list)
    empty_voronoi_chromophore_list = copy.deepcopy(old_chromophore_list)
    for chromo_index, chromo in enumerate(old_chromophore_list):
        empty_cut_off_chromophore_list[chromo_index].neighbours = []
        empty_cut_off_chromophore_list[chromo_index].dissociation_neighbours = []
        empty_voronoi_chromophore_list[chromo_index].neighbours = []
        empty_voronoi_chromophore_list[chromo_index].dissociation_neighbours = []
    sim_dims = [[-axis / 2.0, axis / 2.0] for axis in [AA_morphology_dict[box_length] for box_length
                                                       in ['lx', 'ly', 'lz']]]
    parameter_dict['maximum_hole_hop_distance'] = 10.0
    parameter_dict['maximum_electron_hop_distance'] = 10.0
    old_chromophore_list = obtain_chromophores.determine_neighbours_cut_off(empty_cut_off_chromophore_list,
                                                                            parameter_dict, sim_dims)
    new_chromophore_list = obtain_chromophores.determine_neighbours_voronoi(empty_voronoi_chromophore_list,
                                                                            parameter_dict, sim_dims)
    for list_name in [old_chromophore_list, new_chromophore_list]:
        chromo_ID = 653
        print(chromo_ID)
        print(' '.join(list(map(str, list_name[chromo_ID].AAIDs + [item for sublist in [
            list_name[x[0]].AAIDs for x in list_name[chromo_ID].neighbours] for item in sublist]))))
        print(' '.join(list(map(str, list_name[chromo_ID].AAIDs + [item for sublist in [
            list_name[x[0]].AAIDs for x in list_name[chromo_ID].dissociationNeighbours] for item in sublist]))) + '\n')


def test_write_orca_output(pickle_file):
    # One of the chromophores in the corner is #1198
    R.seed(8585)
    pickle_data = hf.load_pickle(pickle_file)
    AA_morphology_dict = pickle_data[0]
    parameter_dict = pickle_data[3]
    chromophore_list = pickle_data[4]
    for chromo in chromophore_list:
        chromo.neighbours = []
        chromo.dissociation_neighbours = []
    sim_dims = [[-axis / 2.0, axis / 2.0] for axis in [AA_morphology_dict[box_length] for
                                                       box_length in ['lx', 'ly', 'lz']]]
    chromophore_list = obtain_chromophores.determine_neighbours_cut_off(chromophore_list, parameter_dict, sim_dims)
    parameter_dict['output_morph_dir'] = './test_assets/output_files'
    parameter_dict['morphology'] = ''
    execute_ZINDO.create_input_files(chromophore_list, AA_morphology_dict, parameter_dict)


def test_check_periodic_neighbours(pickle_file):
    pickle_data = hf.load_pickle(pickle_file)
    AA_morphology_dict = pickle_data[0]
    parameter_dict = pickle_data[3]
    chromophore_list = pickle_data[4]
    for chromo in chromophore_list:
        chromo.neighbours = []
        chromo.dissociation_neighbours = []
    sim_dims = [[-axis / 2.0, axis / 2.0] for axis in [AA_morphology_dict[box_length] for
                                                       box_length in ['lx', 'ly', 'lz']]]
    chromophore_list = obtain_chromophores.determine_neighbours_voronoi(chromophore_list, parameter_dict, sim_dims)
    chromo_ID = R.randint(0, len(chromophore_list))
    print(chromophore_list[chromo_ID].neighbours)
    print("\nOriginal =", ' '.join(map(str, chromophore_list[chromo_ID].AAIDs)))
    neighbour1_string = "In-image neighbours = "
    neighbour2_string = "Out-of-image neighbours = "
    for [neighbour_ID, image] in chromophore_list[chromo_ID].neighbours:
        if np.array_equal(image, [0, 0, 0]):
            neighbour1_string += ' '.join(map(str, chromophore_list[neighbour_ID].AAIDs)) + ' '
        else:
            neighbour2_string += ' '.join(map(str, chromophore_list[neighbour_ID].AAIDs)) + ' '
    print("\n")
    print(neighbour1_string)
    print("\n")
    print(neighbour2_string)


if __name__ == "__main__":
    # pickle_file = 'test_assets/bilayerBCC/code/bilayerBCC.pickle'
    pickle_file = 'test_assets/p3ht/code/p1-l15-f0.0-p0.1-t1.5-e0.5.pickle'
    # testFindNeighbours(pickleFile)
    test_write_orca_output(pickle_file)
    # testPeriodicNeighbours(pickleFile)
