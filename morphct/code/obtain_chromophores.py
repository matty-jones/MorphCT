import copy
import itertools
import os
import sys
import numpy as np
from collections import defaultdict
from scipy.spatial import Delaunay
from morphct.code import helper_functions as hf


class chromophore:
    def __init__(
        self,
        chromo_ID,
        chromophore_CG_sites,
        CG_morphology_dict,
        AA_morphology_dict,
        CG_to_AAID_master,
        parameter_dict,
        sim_dims,
    ):
        self.ID = chromo_ID
        self.orca_input = "/chromophores/input_orca/single/{:05d}.inp".format(self.ID)
        self.orca_output = "/chromophores/output_orca/single/{:05d}.out".format(self.ID)
        self.CGIDs = chromophore_CG_sites
        # Determine whether this chromophore is a donor or an acceptor, as well
        # as the site types that have been defined as the electronically active
        # in the chromophore
        if CG_morphology_dict is not None:
            # Normal operation
            self.CG_types = sorted(
                list(set([CG_morphology_dict["type"][CGID] for CGID in self.CGIDs]))
            )
            active_CG_sites, self.sub_species = self.obtain_electronic_species(
                chromophore_CG_sites,
                CG_morphology_dict["type"],
                parameter_dict["CG_site_species"],
            )
            self.species = parameter_dict["chromophore_species"][self.sub_species][
                "species"
            ]
            self.reorganisation_energy = parameter_dict["chromophore_species"][
                self.sub_species
            ]["reorganisation_energy"]
            self.VRH_delocalisation = parameter_dict["chromophore_species"][
                self.sub_species
            ]["VRH_delocalisation"]
            # CG_to_AAID_master is a list of dictionaries where each list
            # element corresponds to a new molecule. Firstly, flatten this out
            # so that it becomes a single CG:AAID dictionary
            flattened_CG_to_AAID_master = {
                dict_key: dict_val[1]
                for dictionary in CG_to_AAID_master
                for dict_key, dict_val in dictionary.items()
            }
            # Now, using chromophore_CG_sites as the keys, build up a list of
            # all of the AAIDs in the chromophore, where each element
            # corresponds to each CG site, and then flatten it.
            self.AAIDs = [
                AAID
                for AAIDs in [
                    flattened_CG_to_AAID_master[CGID] for CGID in chromophore_CG_sites
                ]
                for AAID in AAIDs
            ]
            # By using active_CG_sites, determine the AAIDs for
            # the electrically active proportion of the chromophore, so that we
            # can calculate its proper position. Again each element corresponds
            # to each CG site so the list needs to be flattened afterwards.
            electronically_active_AAIDs = [
                AAID
                for AAIDs in [
                    flattened_CG_to_AAID_master[CGID] for CGID in active_CG_sites
                ]
                for AAID in AAIDs
            ]
        else:
            # No fine-graining has been performed by MorphCT, so we know that
            # the input morphology is already atomistic.
            if len(parameter_dict["CG_site_species"]) == 1:
                # If the morphology contains only a single type of electronic
                # species, then the parameter_dict['CG_site_species'] should
                # only have one entry, and we can set all chromophores to be
                # this species.
                active_CG_sites = chromophore_CG_sites
                electronically_active_AAIDs = chromophore_CG_sites
                self.sub_species = list(parameter_dict["CG_site_species"].values())[0]
                self.species = parameter_dict["chromophore_species"][self.sub_species][
                    "species"
                ]
                self.reorganisation_energy = parameter_dict["chromophore_species"][
                    self.sub_species
                ]["reorganisation_energy"]
                self.VRH_delocalisation = parameter_dict["chromophore_species"][
                    self.sub_species
                ]["VRH_delocalisation"]
            elif (len(parameter_dict["CG_site_species"]) == 0) and (
                len(parameter_dict["AA_rigid_body_species"]) > 0
            ):
                # If the CG_site_species have not been specified, then look to
                # the AA_rigid_body_species dictionary to determine which rigid
                # bodies are donors and which are acceptors
                electronically_active_AAIDs = []
                for AAID in chromophore_CG_sites:
                    if AA_morphology_dict["body"][AAID] != -1:
                        electronically_active_AAIDs.append(AAID)
                active_CG_sites = copy.deepcopy(electronically_active_AAIDs)
                # Now work out what the species is:
                for sub_species, rigid_bodies in parameter_dict[
                    "AA_rigid_body_species"
                ].items():
                    if AA_morphology_dict["body"][active_CG_sites[0]] in rigid_bodies:
                        self.sub_species = sub_species
                        self.species = parameter_dict["chromophore_species"][
                            self.sub_species
                        ]["species"]
                        self.reorganisation_energy = parameter_dict[
                            "chromophore_species"
                        ][self.sub_species]["reorganisation_energy"]
                        self.VRH_delocalisation = parameter_dict["chromophore_species"][
                            self.sub_species
                        ]["VRH_delocalisation"]
                        break
                try:
                    self.species
                except AttributeError:
                    for key, val in self.__dict__:
                        print(key, val)
                    raise SystemError(
                        "Chromophore {:d} has no species! Exiting...".format(self.ID)
                    )
            else:
                raise SystemError(
                    "Multiple electronic species defined, but no way to map them"
                    " without a coarse-grained morphology (no CG morph has been given)"
                )
            self.AAIDs = chromophore_CG_sites
        # The position of the chromophore can be calculated easily. Note that
        # here, the `self.image' is the periodic image that the
        # unwrapped_position of the chromophore is located in, relative to the
        # original simulation volume.
        electronically_active_unwrapped_posns = [
            AA_morphology_dict["unwrapped_position"][AAID]
            for AAID in electronically_active_AAIDs
        ]
        electronically_active_types = [
            AA_morphology_dict["type"][AAID] for AAID in electronically_active_AAIDs
        ]
        self.unwrapped_posn, self.posn, self.image = self.obtain_chromophore_COM(
            electronically_active_unwrapped_posns, electronically_active_types, sim_dims
        )
        # A list of the important bonds for this chromophore from the morphology
        # would be useful when determining if a terminating group is already
        # present on this monomer
        self.bonds = self.get_important_bonds(AA_morphology_dict["bond"])
        if CG_morphology_dict is not None:
            # Determine if this chromophore is a repeat unit and therefore will
            # need terminating before orca
            CG_types = set(
                [CG_morphology_dict["type"][CGID] for CGID in chromophore_CG_sites]
            )
            # self.terminate = True if any of the CGTypes in this chromophore
            # are defined as having termination conditions in the parameter file
            self.terminate = any(
                CG_type in CG_types
                for CG_type in [
                    connection[0]
                    for connection in parameter_dict["molecule_terminating_connections"]
                ]
            )
        else:
            try:
                if len(parameter_dict["molecule_terminating_connections"].keys()) == 0:
                    # Small molecules in atomistic morphology therefore no
                    # terminations needed
                    self.terminate = False
            except AttributeError:
                if len(parameter_dict["molecule_terminating_connections"]) == 0:
                    self.terminate = False
            else:
                # No CG morphology, but terminations have been specified, so
                # we're dealing with a polymer
                AA_types = set(
                    [AA_morphology_dict["type"][AAID] for AAID in self.AAIDs]
                )
                self.terminate = any(
                    AA_type in AA_types
                    for AA_type in [
                        connection
                        for connection in parameter_dict[
                            "molecule_terminating_connections"
                        ]
                    ]
                )
        # Now to create a load of placeholder parameters to update later when we
        # have the full list/energy levels.
        # The self.neighbours list contains one element for each chromophore
        # within parameterDict['maximum_hop_distance'] of this one (including
        # periodic boundary conditions). Its format is
        # [[neighbour1_ID, relative_image_of_neighbour1],...]
        self.neighbours = []
        self.dissociation_neighbours = []
        # The molecular orbitals of this chromophore have not yet been
        # calculated, but they will simply be floats.
        self.HOMO = None
        self.HOMO_1 = None
        self.LUMO = None
        self.LUMO_1 = None
        # The neighbour_delta_E and neighbour_TI are lists where each element
        # describes the different in important molecular orbital or transfer
        # integral between this chromophore and each neighbour. The list indices
        # here are the same as in self.neighbours for coherence.
        self.neighbours_delta_E = []
        self.neighbours_TI = []

    def get_important_bonds(self, bond_list):
        important_bonds = []
        for bond in bond_list:
            if (bond[1] in self.AAIDs) and (bond[2] in self.AAIDs):
                important_bonds.append(bond)
        return important_bonds

    def obtain_chromophore_COM(
        self,
        electronically_active_unwrapped_posns,
        electronically_active_types,
        sim_dims,
    ):
        # Calculate the chromophore's position in the morphology (CoM of all
        # atoms in self.AAIDs from AA_morphology_dict)
        chromo_unwrapped_posn = hf.calc_COM(
            electronically_active_unwrapped_posns,
            list_of_atom_types=electronically_active_types,
        )
        chromo_wrapped_posn = copy.deepcopy(chromo_unwrapped_posn)
        chromo_wrapped_image = [0, 0, 0]
        # Now calculate the wrapped position of the chromophore and its image
        for axis in range(3):
            sim_extent = sim_dims[axis][1] - sim_dims[axis][0]
            while chromo_wrapped_posn[axis] < sim_dims[axis][0]:
                chromo_wrapped_posn[axis] += sim_extent
                chromo_wrapped_image[axis] -= 1
            while chromo_wrapped_posn[axis] > sim_dims[axis][1]:
                chromo_wrapped_posn[axis] -= sim_extent
                chromo_wrapped_image[axis] += 1
        return chromo_unwrapped_posn, chromo_wrapped_posn, chromo_wrapped_image

    def obtain_electronic_species(
        self, chromophore_CG_sites, CG_site_types, CG_to_species
    ):
        electronically_active_sites = []
        current_chromophore_species = None
        for CG_site_ID in chromophore_CG_sites:
            site_type = CG_site_types[CG_site_ID]
            site_species = CG_to_species[site_type]
            if site_species.lower() != "none":
                if (current_chromophore_species is not None) and (
                    current_chromophore_species != site_species
                ):
                    raise SystemError(
                        "Problem - multiple electronic species defined in the same "
                        " chromophore. Please modify the chromophore generation code "
                        " to fix this issue for your molecule!"
                    )
                else:
                    current_chromophore_species = site_species
                    electronically_active_sites.append(CG_site_ID)
        return electronically_active_sites, current_chromophore_species

    def get_MO_energy(self):
        if self.species.lower() == "acceptor":
            return self.LUMO
        elif self.species.lower() == "donor":
            return self.HOMO
        else:
            raise Exception("Chromo MUST be donor OR acceptor")


def calculate_chromophores(
    CG_morphology_dict, AA_morphology_dict, CG_to_AAID_master, parameter_dict, sim_dims
):
    # We make the assumption that a chromophore consists of one of each of the
    # CG site types described by the same template file. For instance, if we
    # have 3 sites 'A', 'B' and 'C' described in one file and one site 'D'
    # described in another file then there are two chromophores species
    # described by A-B-C and D. This will be treated automatically because the
    # D's shouldn't be bonded to anything in the CGMorphologyDict if they are
    # small molecules.
    # Therefore, we need to assign each CG site in the morphology to a
    # particular chromophore, so first, it's important to generate a
    # `neighbour_list' of all bonded atoms
    print("Determining chromophores in the system...")
    bonded_atoms = hf.obtain_bonded_list(CG_morphology_dict["bond"])
    chromophore_list = [i for i in range(len(CG_morphology_dict["type"]))]
    for CG_site_ID, chromophore_ID in enumerate(chromophore_list):
        CG_site_type = CG_morphology_dict["type"][CG_site_ID]
        types_in_this_chromophore = [CG_site_type]
        chromophore_list, types_in_this_chromophore = update_chromophores(
            CG_site_ID,
            chromophore_list,
            bonded_atoms,
            CG_morphology_dict["type"],
            types_in_this_chromophore,
            parameter_dict,
        )
    chromophore_data = {}
    for atom_ID, chromo_ID in enumerate(chromophore_list):
        if chromo_ID not in list(chromophore_data.keys()):
            chromophore_data[chromo_ID] = [atom_ID]
        else:
            chromophore_data[chromo_ID].append(atom_ID)
    # Now rename the chromophore IDs so that they increment sensibly (they will
    # be used later for the orca files)
    old_keys = sorted(chromophore_data.keys())
    for new_key, old_key in enumerate(old_keys):
        chromophore_data[new_key] = chromophore_data.pop(old_key)
    print(
        "{:d} chromophores successfully identified!".format(
            len(list(chromophore_data.keys()))
        )
    )
    # Now let's create a list of all the chromophore instances which contain all
    # of the information we could ever want about them.
    chromophore_instances = []
    for chromo_ID, chromophore_CG_sites in chromophore_data.items():
        print(
            "\rCalculating properties of chromophore {:05d} of {:05d}...".format(
                chromo_ID, len(list(chromophore_data.keys())) - 1
            ),
            end=" ",
        )
        if sys.stdout is not None:
            sys.stdout.flush()
        chromophore_instances.append(
            chromophore(
                chromo_ID,
                chromophore_CG_sites,
                CG_morphology_dict,
                AA_morphology_dict,
                CG_to_AAID_master,
                parameter_dict,
                sim_dims,
            )
        )
    print("")
    return chromophore_instances


def calculate_chromophores_AA(
    CG_morphology_dict,
    AA_morphology_dict,
    CG_to_AAID_master,
    parameter_dict,
    sim_dims,
    rigid_bodies=None,
):
    # If rigid_bodies == None:
    # This function works in the same way as the coarse-grained version above,
    # except this one iterates through the AA bonds instead. This is FAR SLOWER
    # and so shouldn't be done, except in the case where the coarse-grained
    # morphology does not exist (because we started with an atomistic morphology
    # and are only interested in running KMC on it)
    # If rigid_bodies == AA_morphology_dict['body']:
    # This function uses the rigid bodies specified in
    # parameter_dict['AA_rigid_body_species'], and those which have not been
    # specified by iterating through the AA bond list, to determine the
    # chromophores in the system. This is the slowest way to calculate
    # chromophores, but is useful for systems such as BDT-TPD, where there are
    # multiple chromophores of differing species present in the same molecule.
    # As above, this code will only run if an atomistic morphology has been
    # input to MorphCT. If it is coarse-grained, the CG-based
    # "calculate_chromophore" function will be used, and will also be a lot
    # faster.
    # The parameter_dict['AA_rigid_body_species'] is a dictionary with two keys,
    # 'donor' or 'acceptor'. Each element in the value list corresponds to a new
    # chromophore. These aren't the only atoms that belong to this chromophore,
    # however - there might be a bunch of aliphatic/flexible atoms that are
    # connected, so we need to make sure that we add those too.
    print("Determining chromophores in the system...")
    bonded_atoms = hf.obtain_bonded_list(AA_morphology_dict["bond"])
    chromophore_list = [i for i in range(len(AA_morphology_dict["type"]))]
    for AA_site_ID, chromophore_ID in enumerate(chromophore_list):
        AA_site_type = AA_morphology_dict["type"][AA_site_ID]
        chromophore_list = update_chromophores_AA(
            AA_site_ID, chromophore_list, bonded_atoms, parameter_dict, rigid_bodies
        )
    chromophore_data = {}
    for atom_ID, chromo_ID in enumerate(chromophore_list):
        if chromo_ID not in list(chromophore_data.keys()):
            chromophore_data[chromo_ID] = [atom_ID]
        else:
            chromophore_data[chromo_ID].append(atom_ID)
    # Now rename the chromophore IDs so that they increment sensibly (they will
    # be used later for the orca files)
    old_keys = sorted(chromophore_data.keys())
    for new_key, old_key in enumerate(old_keys):
        chromophore_data[new_key] = chromophore_data.pop(old_key)
    print(
        "{:d} chromophores successfully identified!".format(
            len(list(chromophore_data.keys()))
        )
    )
    # Now let's create a list of all the chromophore instances which contain all
    # of the information we could ever want about them.
    chromophore_instances = []
    for chromo_ID, chromophore_CG_sites in chromophore_data.items():
        print(
            "\rCalculating properties of chromophore {:05d} of {:05d}...".format(
                chromo_ID, len(list(chromophore_data.keys())) - 1
            ),
            end=" ",
        )
        if sys.stdout is not None:
            sys.stdout.flush()
        chromophore_instances.append(
            chromophore(
                chromo_ID,
                chromophore_CG_sites,
                CG_morphology_dict,
                AA_morphology_dict,
                CG_to_AAID_master,
                parameter_dict,
                sim_dims,
            )
        )
    print("")
    return chromophore_instances


def update_chromophores(
    atom_ID,
    chromophore_list,
    bonded_atoms,
    CG_type_list,
    types_in_this_chromophore,
    parameter_dict,
):
    # Recursively add all neighbours of atom number atom_ID to this chromophore,
    # providing the same type does not already exist in it
    try:
        for bonded_atom in bonded_atoms[atom_ID]:
            bonded_type = CG_type_list[bonded_atom]
            # First check that the bonded_atom's type is not already in this
            # chromophore.
            # Also, check that the type to be added is of the same electronic
            # species as the ones added previously, or == 'None'
            if (bonded_type not in types_in_this_chromophore) and (
                (parameter_dict["CG_site_species"][bonded_type].lower() == "none")
                or (
                    parameter_dict["CG_site_species"][bonded_type].lower()
                    == list(
                        set(
                            [
                                parameter_dict["CG_site_species"][x].lower()
                                for x in types_in_this_chromophore
                            ]
                        )
                    )[0]
                )
            ):
                # If the atomID of the bonded atom is larger than that of the
                # current one, update the bonded atom's ID to the current one's
                # to put it in this chromophore, then iterate through all of the
                # bonded atom's neighbours
                if chromophore_list[bonded_atom] > chromophore_list[atom_ID]:
                    chromophore_list[bonded_atom] = chromophore_list[atom_ID]
                    types_in_this_chromophore.append(bonded_type)
                    chromophore_list, types_in_this_chromophore = update_chromophores(
                        bonded_atom,
                        chromophore_list,
                        bonded_atoms,
                        CG_type_list,
                        types_in_this_chromophore,
                        parameter_dict,
                    )
                # If the atomID of the current atom is larger than that of the
                # bonded one, update the current atom's ID to the bonded one's
                # to put it in this chromophore, then iterate through all of the
                # current atom's neighbours
                elif chromophore_list[bonded_atom] < chromophore_list[atom_ID]:
                    chromophore_list[atom_ID] = chromophore_list[bonded_atom]
                    types_in_this_chromophore.append(CG_type_list[atom_ID])
                    chromophore_list, types_in_this_chromophore = update_chromophores(
                        atom_ID,
                        chromophore_list,
                        bonded_atoms,
                        CG_type_list,
                        types_in_this_chromophore,
                        parameter_dict,
                    )
                # Else: both the current and the bonded atom are already known
                # to be in this chromophore, so we don't have to do anything
                # else.
    except KeyError:
        # This means that there are no bonded CG sites (i.e. it's a single
        # chromophore)
        pass
    return chromophore_list, types_in_this_chromophore


def update_chromophores_AA(
    atom_ID, chromophore_list, bonded_atoms, parameter_dict, rigid_bodies=None
):
    # This version of the update chromophores function does not check for CG
    # site types, instead just adding all bonded atoms. Therefore it should only
    # be used in the case of already-atomistic morphologies (no CG morph
    # specified) containing ONLY small molecules
    try:
        for bonded_atom in bonded_atoms[atom_ID]:
            if rigid_bodies is not None:
                # Skip if the bonded atom belongs to a different rigid body
                if (
                    (rigid_bodies[bonded_atom] != -1) and (rigid_bodies[atom_ID] != -1)
                ) and (rigid_bodies[bonded_atom] != rigid_bodies[atom_ID]):
                    continue
            # If the atomID of the bonded atom is larger than that of the
            # current one, update the bonded atom's ID to the current one's to
            # put it in this chromophore, then iterate through all of the bonded
            # atom's neighbours
            if chromophore_list[bonded_atom] > chromophore_list[atom_ID]:
                chromophore_list[bonded_atom] = chromophore_list[atom_ID]
                chromophore_list = update_chromophores_AA(
                    bonded_atom,
                    chromophore_list,
                    bonded_atoms,
                    parameter_dict,
                    rigid_bodies,
                )
            # If the atomID of the current atom is larger than that of the
            # bonded one, update the current atom's ID to the bonded one's to
            # put it in this chromophore, then iterate through all of the
            # current atom's neighbours
            elif chromophore_list[bonded_atom] < chromophore_list[atom_ID]:
                chromophore_list[atom_ID] = chromophore_list[bonded_atom]
                chromophore_list = update_chromophores_AA(
                    atom_ID,
                    chromophore_list,
                    bonded_atoms,
                    parameter_dict,
                    rigid_bodies,
                )
            # Else: both the current and the bonded atom are already known to be
            # in this chromophore, so we don't have to do anything else.
    except KeyError:
        # This means that there are no bonded CG sites (i.e. it's a single
        # chromophore)
        pass
    return chromophore_list


def create_super_cell(chromophore_list, box_size):
    for chromophore in chromophore_list:
        chromophore.super_cell_posns = []
        chromophore.super_cell_images = []
        for x_image in range(-1, 2):
            for y_image in range(-1, 2):
                for z_image in range(-1, 2):
                    chromophore.super_cell_posns.append(
                        np.array(chromophore.posn)
                        + (np.array([x_image, y_image, z_image]) * (np.array(box_size)))
                    )
                    chromophore.super_cell_images.append(
                        np.array([x_image, y_image, z_image])
                    )
    return chromophore_list


def get_voronoi_neighbours(tri, chromo_list):
    n_list = defaultdict(set)
    for p in tri.vertices:
        for i, j in itertools.permutations(p, 2):
            n_list[chromo_list[i].periodic_ID].add(chromo_list[j].periodic_ID)
    return n_list


class super_cell_chromo:
    def __init__(self):
        self.species = None
        self.original_ID = None
        self.periodic_ID = None
        self.position = None
        self.image = None


def update_chromophore_list_voronoi(
    IDs_to_update, super_cell_chromos, neighbour_IDs, chromophore_list, sim_dims
):
    # IDs to Update is a list of the periodic chromophores with the image
    # [0, 0, 0]
    for periodic_ID in IDs_to_update:
        # Obtain the real chromophore corresponding to this periodic_ID
        chromophore1 = chromophore_list[super_cell_chromos[periodic_ID].original_ID]
        assert np.array_equal(super_cell_chromos[periodic_ID].image, [0, 0, 0])
        # Get latest neighbour information
        chromo1neighbour_IDs = [
            neighbour_data[0] for neighbour_data in chromophore1.neighbours
        ]
        chromo1dissociation_neighbour_IDs = [
            neighbour_data[0] for neighbour_data in chromophore1.dissociation_neighbours
        ]
        for neighbour_periodic_ID in neighbour_IDs[periodic_ID]:
            neighbour_super_cell_chromo = super_cell_chromos[neighbour_periodic_ID]
            chromophore2 = chromophore_list[neighbour_super_cell_chromo.original_ID]
            chromo2neighbour_IDs = [
                neighbour_data[0] for neighbour_data in chromophore2.neighbours
            ]
            chromo2dissociation_neighbour_IDs = [
                neighbour_data[0]
                for neighbour_data in chromophore2.dissociation_neighbours
            ]
            relative_image = neighbour_super_cell_chromo.image
            if chromophore1.species == chromophore2.species:
                if chromophore2.ID not in chromo1neighbour_IDs:
                    chromophore1.neighbours.append(
                        [chromophore2.ID, list(np.array(relative_image))]
                    )
                    chromophore1.neighbours_delta_E.append(None)
                    chromophore1.neighbours_TI.append(None)
                    chromo1neighbour_IDs.append(chromophore2.ID)
                if chromophore1.ID not in chromo2neighbour_IDs:
                    chromophore2.neighbours.append(
                        [chromophore1.ID, list(-np.array(relative_image))]
                    )
                    chromophore2.neighbours_delta_E.append(None)
                    chromophore2.neighbours_TI.append(None)
                    chromo2neighbour_IDs.append(chromophore1.ID)
            else:
                if chromophore2.ID not in chromo1dissociation_neighbour_IDs:
                    chromophore1.dissociation_neighbours.append(
                        [chromophore2.ID, list(np.array(relative_image))]
                    )
                    chromo1dissociation_neighbour_IDs.append(chromophore2.ID)
                if chromophore1.ID not in chromo2dissociation_neighbour_IDs:
                    chromophore2.dissociation_neighbours.append(
                        [chromophore1.ID, list(-np.array(relative_image))]
                    )
                    chromo2dissociation_neighbour_IDs.append(chromophore1.ID)
    return chromophore_list


def determine_neighbours_voronoi(chromophore_list, parameter_dict, sim_dims):
    box_size = [axis[1] - axis[0] for axis in sim_dims]
    # First create the supercell
    super_cell = create_super_cell(chromophore_list, box_size)
    donor_chromos = []
    acceptor_chromos = []
    all_chromos = []
    chromo_index = 0
    for chromophore in super_cell:
        for index, position in enumerate(chromophore.super_cell_posns):
            chromo = super_cell_chromo()
            chromo.species = chromophore.species
            chromo.original_ID = chromophore.ID
            chromo.periodic_ID = chromo_index
            chromo.position = position
            chromo.image = chromophore.super_cell_images[index]
            chromo_index += 1
            if chromophore.species.lower() == "donor":
                donor_chromos.append(chromo)
            elif chromophore.species.lower() == "acceptor":
                acceptor_chromos.append(chromo)
            all_chromos.append(chromo)
    # Now obtain the positions and send them to the Delaunay Triangulation
    # Then get the voronoi neighbours
    all_positions = [chromo.position for chromo in all_chromos]
    # Initialise the neighbour dictionaries
    all_neighbours = {}
    # Update the relevant neighbour dictionaries if we have the right
    # chromophore types in the system. Also log the chromophoreIDs from the
    # original simulation volume (non-periodic). Chromophores in the original
    # simulation volume will be every 27th (there are 27 periodic images in the
    # triple range(-1,2)), beginning from #13 ((0, 0, 0) is the thirteenth
    # element of the triple range(-1,2)) up to the length of the list in
    # question.
    original_all_chromo_IDs = []
    try:
        if parameter_dict["permit_hops_through_opposing_chromophores"]:
            # Need to only consider the neighbours of like chromophore species
            donor_positions = [chromo.position for chromo in donor_chromos]
            acceptor_positions = [chromo.position for chromo in acceptor_chromos]
            donor_neighbours = {}
            acceptor_neighbours = {}
            original_donor_chromo_IDs = []
            original_acceptor_chromo_IDs = []
            for chromophore in all_chromos:
                if np.array_equal(chromophore.image, [0, 0, 0]):
                    original_all_chromo_IDs.append(chromophore.periodic_ID)
                    if chromophore.species.lower() == "donor":
                        original_donor_chromo_IDs.append(chromophore.periodic_ID)
                    elif chromophore.species.lower() == "acceptor":
                        original_acceptor_chromo_IDs.append(chromophore.periodic_ID)
            if len(donor_positions) > 0:
                print("Calculating Neighbours of donor Moieties")
                donor_neighbours = get_voronoi_neighbours(
                    Delaunay(donor_positions), donor_chromos
                )
                print("Updating the chromophore list for donor chromos")
                chromophore_list = update_chromophore_list_voronoi(
                    original_donor_chromo_IDs,
                    all_chromos,
                    donor_neighbours,
                    chromophore_list,
                    sim_dims,
                )
            if len(acceptor_positions) > 0:
                print("Calculating Neighbours of acceptor Moieties")
                acceptor_neighbours = get_voronoi_neighbours(
                    Delaunay(acceptor_positions), acceptor_chromos
                )
                print("Updating the chromophore list for acceptor chromos")
                chromophore_list = update_chromophore_list_voronoi(
                    original_acceptor_chromo_IDs,
                    all_chromos,
                    acceptor_neighbours,
                    chromophore_list,
                    sim_dims,
                )
        else:
            raise KeyError
    except KeyError:
        # Default behaviour - carriers are blocked by the opposing species
        for chromophore in all_chromos:
            if np.array_equal(chromophore.image, [0, 0, 0]):
                original_all_chromo_IDs.append(chromophore.periodic_ID)
    print("Calculating Neighbours of All Moieties")
    all_neighbours = get_voronoi_neighbours(Delaunay(all_positions), all_chromos)
    print("Updating the chromophore list for dissociation neighbours")
    chromophore_list = update_chromophore_list_voronoi(
        original_all_chromo_IDs, all_chromos, all_neighbours, chromophore_list, sim_dims
    )
    return chromophore_list


def determine_neighbours_cut_off(chromophore_list, parameter_dict, sim_dims):
    for chromophore1 in chromophore_list:
        print(
            "\rIdentifying neighbours of chromophore {:05d} of {:05d}...".format(
                chromophore1.ID, len(chromophore_list) - 1
            ),
            end=" ",
        )
        if sys.stdout is not None:
            sys.stdout.flush()
        for chromophore2 in chromophore_list:
            # Skip if chromo2 is chromo1
            if chromophore1.ID == chromophore2.ID:
                continue
            delta_posn = chromophore2.posn - chromophore1.posn
            relative_image_of_chromo2 = [0, 0, 0]
            # Consider periodic boundary conditions
            for axis in range(3):
                half_box_length = (sim_dims[axis][1] - sim_dims[axis][0]) / 2.0
                while delta_posn[axis] > half_box_length:
                    delta_posn[axis] -= sim_dims[axis][1] - sim_dims[axis][0]
                    relative_image_of_chromo2[axis] -= 1
                while delta_posn[axis] < -half_box_length:
                    delta_posn[axis] += sim_dims[axis][1] - sim_dims[axis][0]
                    relative_image_of_chromo2[axis] += 1
            separation = np.linalg.norm(delta_posn)
            # If proximity is within tolerance, add these chromophores as
            # neighbours. Base check is against the maximum of the donor and
            # acceptor hop distances. A further separation check is made if the
            # chromophores are the same type to make sure we don't exceed the
            # maximum specified hop distance for the carrier type.
            if separation <= max(
                [
                    parameter_dict["maximum_hole_hop_distance"],
                    parameter_dict["maximum_electron_hop_distance"],
                ]
            ):
                # Only add the neighbours if they haven't already been added so
                # far
                chromo1neighbour_IDs = [
                    neighbour_data[0] for neighbour_data in chromophore1.neighbours
                ]
                chromo2neighbour_IDs = [
                    neighbour_data[0] for neighbour_data in chromophore2.neighbours
                ]
                chromo1dissociation_neighbour_IDs = [
                    neighbour_data[0]
                    for neighbour_data in chromophore1.dissociation_neighbours
                ]
                chromo2dissociation_neighbour_IDs = [
                    neighbour_data[0]
                    for neighbour_data in chromophore2.dissociation_neighbours
                ]
                # Also, make the delta_E and the T_ij lists as long as the
                # neighbour lists for easy access later
                if chromophore1.species == chromophore2.species:
                    if (
                        (chromophore1.species.lower() == "donor")
                        and (separation >= parameter_dict["maximum_hole_hop_distance"])
                    ) or (
                        (chromophore1.species.lower() == "acceptor")
                        and (
                            separation
                            >= parameter_dict["maximum_electron_hop_distance"]
                        )
                    ):
                        continue
                    if chromophore2.ID not in chromo1neighbour_IDs:
                        chromophore1.neighbours.append(
                            [chromophore2.ID, relative_image_of_chromo2]
                        )
                        chromophore1.neighbours_delta_E.append(None)
                        chromophore1.neighbours_TI.append(None)
                    if chromophore1.ID not in chromo2neighbour_IDs:
                        chromophore2.neighbours.append(
                            [
                                chromophore1.ID,
                                list(-np.array(relative_image_of_chromo2)),
                            ]
                        )
                        chromophore2.neighbours_delta_E.append(None)
                        chromophore2.neighbours_TI.append(None)
                else:
                    # NOTE: Modifying this so that only dissociation neigbours in the
                    # same periodic image are considered.
                    if (chromophore2.ID not in chromo2dissociation_neighbour_IDs) and (
                        np.all(np.isclose(relative_image_of_chromo2, [0, 0, 0]))
                    ):
                        chromophore1.dissociation_neighbours.append(
                            [chromophore2.ID, [0, 0, 0]]
                        )
                    if (chromophore1.ID not in chromo1dissociation_neighbour_IDs) and (
                        np.all(np.isclose(relative_image_of_chromo2, [0, 0, 0]))
                    ):
                        chromophore2.dissociation_neighbours.append(
                            [chromophore1.ID, [0, 0, 0]]
                        )
    print("")
    return chromophore_list


def chromo_sort(chromophore_list):
    for index, chromo in enumerate(chromophore_list):
        if index != chromo.ID:
            print(
                "Inconsistency found in the ordering of the chromophore_list, rewriting"
                " the chromophore_list in the correct order..."
            )
            new_chromophore_list = []
            for chromo in chromophore_list:
                new_chromophore_list.append(0)
            for chromo in chromophore_list:
                new_chromophore_list[chromo.ID] = chromo
            chromophore_list = new_chromophore_list
            return chromophore_list
    return chromophore_list


def main(
    AA_morphology_dict,
    CG_morphology_dict,
    CG_to_AAID_master,
    parameter_dict,
    chromophore_list,
):
    sim_dims = [
        [-AA_morphology_dict["lx"] / 2.0, AA_morphology_dict["lx"] / 2.0],
        [-AA_morphology_dict["ly"] / 2.0, AA_morphology_dict["ly"] / 2.0],
        [-AA_morphology_dict["lz"] / 2.0, AA_morphology_dict["lz"] / 2.0],
    ]
    if len(parameter_dict["CG_to_template_dirs"]) > 0:
        # Normal operation using the coarse-grained morphology
        chromophore_list = calculate_chromophores(
            CG_morphology_dict,
            AA_morphology_dict,
            CG_to_AAID_master,
            parameter_dict,
            sim_dims,
        )
    elif (len(parameter_dict["CG_site_species"]) == 1) and (
        len(parameter_dict["AA_rigid_body_species"]) == 0
    ):
        # Small molecule system with only one electronic species
        chromophore_list = calculate_chromophores_AA(
            CG_morphology_dict,
            AA_morphology_dict,
            CG_to_AAID_master,
            parameter_dict,
            sim_dims,
        )
    else:
        # Other system, with electronically active species specified as rigid
        # bodies using AA_rigid_body_species in parameter file
        chromophore_list = calculate_chromophores_AA(
            CG_morphology_dict,
            AA_morphology_dict,
            CG_to_AAID_master,
            parameter_dict,
            sim_dims,
            rigid_bodies=AA_morphology_dict["body"],
        )
    chromophore_list = chromo_sort(chromophore_list)
    if parameter_dict["use_voronoi_neighbours"] is True:
        chromophore_list = determine_neighbours_voronoi(
            chromophore_list, parameter_dict, sim_dims
        )
    else:
        chromophore_list = determine_neighbours_cut_off(
            chromophore_list, parameter_dict, sim_dims
        )
    # Now we have updated the chromophore_list, rewrite the pickle with this new
    # information.
    pickle_name = os.path.join(
        parameter_dict["output_morphology_directory"],
        "code",
        "".join([os.path.splitext(parameter_dict["morphology"])[0], ".pickle"]),
    )
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
