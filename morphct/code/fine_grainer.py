import copy
import os
import sys
import numpy as np
from morphct.code import helper_functions as hf


class morphology:
    def __init__(
        self, morphology_xml, morphology_name, parameter_dict, chromophore_list
    ):
        # Need to save the parameter_dict in full as well as its component
        # values because we're going to update the parameterDict with the new
        # type mappings by the end of this code module.
        self.parameter_dict = parameter_dict
        # Import parameters from the parXX.py
        for key, value in parameter_dict.items():
            self.__dict__[key] = value
        self.xml_path = morphology_xml
        self.morphology_name = morphology_name
        # self.inputSigma is the `compression value' in Angstroms that has been
        # used to scale the morphology
        # E.G. the P3HT template uses sigma = 1, but the Marsh morphologies use
        # sigma = 3.
        self.CG_dictionary = hf.load_morphology_xml(
            self.xml_path, sigma=self.input_sigma
        )
        self.CG_dictionary = hf.add_unwrapped_positions(self.CG_dictionary)
        self.chromophore_list = chromophore_list

    def analyse_morphology(self):
        # Split the morphology into individual molecules
        print("Finding molecules...")
        molecule_IDs, molecule_lengths = self.split_molecules()
        rolling_AA_index = 0
        CG_morphology_dict = {}
        AA_morphology_dict = {}
        # Set the AAMorphology and CGMorphology system sizes to the same as the
        # input file system size
        for box_dimension in ["lx", "ly", "lz", "xy", "xz", "yz"]:
            CG_morphology_dict[box_dimension] = self.CG_dictionary[box_dimension]
            AA_morphology_dict[box_dimension] = self.CG_dictionary[box_dimension]
        CG_to_AAID_master = []  # This is a list of dictionaries. Elements in
        # the list correspond to molecules (useful for splitting out individual
        # molecules for the xyz conversion) within the element, the dictionary
        # key is the CG site, the value is a list containing the CG type (e.g.
        # 'thio') as the first element and then another list of all the AAIDs
        # corresponding to that CG site as the second element.

        # If no CG_to_template info is present in the parameter dict, then we
        # can assume that the morphology is already fine-grained and so we can
        # just return the important information and skip this module
        if len(self.CG_to_template_dirs) == 0:
            print(
                "No CG to AA data found in parameter file - the morphology is already"
                " fine-grained! Skipping this module..."
            )
            # Write the xml file and create the pickle
            print("Writing xml file...")
            AA_file_name = os.path.join(
                self.output_morph_dir,
                self.morphology_name,
                "morphology",
                "".join([self.morphology_name, ".xml"]),
            )
            atomistic_morphology = hf.add_unwrapped_positions(self.CG_dictionary)
            # Now write the morphology xml
            hf.write_morphology_xml(atomistic_morphology, AA_file_name)
            # And finally write the pickle with the CGDictionary as None (to
            # indicate to MorphCT that no fine-graining has taken place), but
            # the other parameters assigned as required.
            pickle_location = os.path.join(
                self.output_morph_dir,
                self.morphology_name,
                "code",
                "".join([self.morphology_name, ".pickle"]),
            )
            hf.write_pickle(
                (
                    atomistic_morphology,
                    None,
                    None,
                    self.parameter_dict,
                    self.chromophore_list,
                ),
                pickle_location,
            )
            return (
                atomistic_morphology,
                None,
                None,
                self.parameter_dict,
                self.chromophore_list,
            )

        # Create a ghost particle dictionary to be added at the end of the
        # morphology. This way, we don't mess up the order of atoms later on
        # when trying to split back up into individual molecules and monomers.
        # The ghost dictionary contains all of the type T and type X particles
        # that will anchor the thiophene rings to the CG COM positions.
        ghost_dictionary = {
            "position": [],
            "image": [],
            "unwrapped_position": [],
            "mass": [],
            "diameter": [],
            "type": [],
            "body": [],
            "bond": [],
            "angle": [],
            "dihedral": [],
            "improper": [],
            "charge": [],
        }

        # Need to check for atom-type conflicts and suitably increment the
        # type indices if more than one molecule type is being fine-grained
        new_type_mappings = self.get_new_type_mappings(
            self.CG_to_template_dirs, self.CG_to_template_force_fields
        )
        # Need to update the self.parameterDict, which will be rewritten at the
        # end of this module
        self.parameter_dict["new_type_mappings"] = new_type_mappings
        molecule = []
        unique_mappings = []
        CG_sites, mappings = hf.parallel_sort(
            list(new_type_mappings.keys()), list(new_type_mappings.values())
        )
        for index, mapping in enumerate(mappings):
            if mapping not in unique_mappings:
                molecule.append([])
                unique_mappings.append(mapping)
            molecule[-1].append(CG_sites[index])
        printExplanation = True
        for index, CG_sites in enumerate(molecule):
            printMol = True
            initial_atoms, final_atoms = hf.parallel_sort(
                list(unique_mappings[index].keys()),
                list(unique_mappings[index].values()),
            )
            for index, initial_atom in enumerate(initial_atoms):
                if initial_atom == final_atoms[index]:
                    continue
                if printExplanation is True:
                    print(
                        "The following atom types have been remapped due to conflicting"
                        " typenames in the atomistic templates:"
                    )
                    printExplanation = False
                if printMol is True:
                    print(
                        "Atom types belonging the molecule described by",
                        "".join([repr(CG_sites), ":"]),
                    )
                    printMol = False
                print(initialAtom, "--->", finalAtoms[index])
        print("Adding", len(molecule_IDs), "molecules to the system...")
        for molecule_number in range(len(molecule_IDs)):
            print("Adding molecule number", molecule_number, "\r", end=" ")
            if sys.stdout is not None:
                sys.stdout.flush()
            # Obtain the AA dictionary for each molecule using the
            # fine-graining procedure
            AA_molecule_dict, C_gto_AAIDs, ghost_dictionary = atomistic(
                molecule_number,
                molecule_IDs[molecule_number],
                self.CG_dictionary,
                molecule_lengths,
                rolling_AA_index,
                ghost_dictionary,
                self.parameter_dict,
            ).return_data()
            CG_to_AAID_master.append(C_gto_AAIDs)
            # Update the morphology dictionaries with this new molecule
            for key in list(self.CG_dictionary.keys()):
                if key not in [
                    "lx",
                    "ly",
                    "lz",
                    "xy",
                    "xz",
                    "yz",
                    "time_step",
                    "dimensions",
                ]:
                    if key not in list(AA_morphology_dict.keys()):
                        AA_morphology_dict[key] = AA_molecule_dict[key]
                    else:
                        AA_morphology_dict[key] += AA_molecule_dict[key]
            rolling_AA_index += len(AA_molecule_dict["type"])
        # Now add the ghost dictionary to the end of the morphology file
        # total_number_of_atoms should be == rolling_AA_index, but don't want
        # to take any chances.
        total_number_of_atoms = len(AA_morphology_dict["type"])
        # Add in the wrapped positions of the ghosts. Need to know sim dims for
        # this
        for key in ["lx", "ly", "lz", "xy", "xz", "yz"]:
            ghost_dictionary[key] = AA_morphology_dict[key]
        ghost_dictionary = hf.add_wrapped_positions(ghost_dictionary)
        for key in ["lx", "ly", "lz", "xy", "xz", "yz"]:
            ghost_dictionary.pop(key)
        # The real atoms that the ghost particles are bonded to are already
        # correct and no longer need to be changed.
        # However, the ghost particles themselves will have the wrong indices
        # if we were to add them to the system directly.
        # Therefore, increment all of the ghost bond indices that begin with a
        # * (ghost particle) by the total number of atoms already in the
        # system.
        for bond_no, bond in enumerate(ghost_dictionary["bond"]):
            if str(bond[1])[0] == "*":
                ghost_dictionary["bond"][bond_no][1] = (
                    int(bond[1][1:]) + total_number_of_atoms
                )
            if str(bond[2])[0] == "*":
                ghost_dictionary["bond"][bond_no][2] = (
                    int(bond[2][1:]) + total_number_of_atoms
                )
        # Now append all ghosts to morphology
        for key in list(ghost_dictionary.keys()):
            AA_morphology_dict[key] += ghost_dictionary[key]
        # Finally, update the number of atoms
        AA_morphology_dict["natoms"] += len(ghost_dictionary["type"])
        print("\n")
        # Now write the xml file and create the pickle
        print("Writing xml file...")
        AA_file_name = os.path.join(
            self.output_morph_dir,
            self.morphology_name,
            "morphology",
            "".join([self.morphology_name, ".xml"]),
        )
        # Replace the `positions' with the `unwrapped_positions' ready for
        # writing
        AA_morphology_dict = hf.replace_wrapped_positions(AA_morphology_dict)
        # Update the additional_constraints that we put in by checking all of
        # the constraints have the correct names before writing
        AA_morphology_dict = hf.check_constraint_names(AA_morphology_dict)
        # Now write the morphology xml
        hf.write_morphology_xml(AA_morphology_dict, AA_file_name)
        # And finally write the pickle
        pickle_location = os.path.join(
            self.output_morph_dir,
            self.morphology_name,
            "code",
            "".join([self.morphology_name, ".pickle"]),
        )
        hf.write_pickle(
            (
                AA_morphology_dict,
                self.CG_dictionary,
                CG_to_AAID_master,
                self.parameter_dict,
                self.chromophore_list,
            ),
            pickle_location,
        )
        return (
            AA_morphology_dict,
            self.CG_dictionary,
            CG_to_AAID_master,
            self.parameter_dict,
            self.chromophore_list,
        )

    def split_molecules(self):
        # Split the full morphology into individual molecules
        molecule_AAIDs = []
        molecule_lengths = []
        # Create a lookup table `neighbour list' for all connected atoms called
        # {bonded_atoms}
        bonded_atoms = hf.obtain_bonded_list(self.CG_dictionary["bond"])
        molecule_list = [i for i in range(len(self.CG_dictionary["type"]))]
        # Recursively add all atoms in the neighbour list to this molecule
        for mol_ID in range(len(molecule_list)):
            molecule_list = self.update_molecule(mol_ID, molecule_list, bonded_atoms)
        # Create a dictionary of the molecule data
        molecule_data = {}
        for atom_ID in range(len(self.CG_dictionary["type"])):
            if molecule_list[atom_ID] not in molecule_data:
                molecule_data[molecule_list[atom_ID]] = [atom_ID]
            else:
                molecule_data[molecule_list[atom_ID]].append(atom_ID)
        # Return the list of AAIDs and the lengths of the molecules
        for molecule_ID in list(molecule_data.keys()):
            molecule_AAIDs.append(sorted(molecule_data[molecule_ID]))
            molecule_lengths.append(len(molecule_data[molecule_ID]))
        return molecule_AAIDs, molecule_lengths

    def update_molecule(self, atom_ID, molecule_list, bonded_atoms):
        # Recursively add all neighbours of atom number atomID to this molecule
        try:
            for bonded_atom in bonded_atoms[atom_ID]:
                # If the moleculeID of the bonded atom is larger than that of
                # the current one, update the bonded atom's ID to the current
                # one's to put it in this molecule, then iterate through all of
                # the bonded atom's neighbours
                if molecule_list[bonded_atom] > molecule_list[atom_ID]:
                    molecule_list[bonded_atom] = molecule_list[atom_ID]
                    molecule_list = self.update_molecule(
                        bonded_atom, molecule_list, bonded_atoms
                    )
                # If the moleculeID of the current atom is larger than that of
                # the bonded one, update the current atom's ID to the bonded
                # one's to put it in this molecule, then iterate through all of
                # the current atom's neighbours
                elif molecule_list[bonded_atom] < molecule_list[atom_ID]:
                    molecule_list[atom_ID] = molecule_list[bonded_atom]
                    molecule_list = self.update_molecule(
                        atom_ID, molecule_list, bonded_atoms
                    )
                # Else: both the current and the bonded atom are already known
                # to be in this molecule, so we don't have to do anything else.
        except KeyError:
            # This means that there are no bonded CG sites (i.e. it's a single
            # molecule)
            pass
        return molecule_list

    def get_new_type_mappings(self, CG_to_template_dirs, CG_to_template_force_fields):
        force_field_locations = []
        force_field_mappings = []
        morphology_atom_types = []
        CG_to_template_mappings = {}
        for CG_site, directory in CG_to_template_dirs.items():
            FF_loc = os.path.join(directory, CG_to_template_force_fields[CG_site])
            if FF_loc not in force_field_locations:
                force_field_locations.append(FF_loc)
        for FF_loc in force_field_locations:
            mapping_for_this_FF = {}
            force_field = hf.load_FF_xml(FF_loc)
            for lj_interaction in force_field["lj"]:
                atom_type = lj_interaction[0]
                while atom_type in morphology_atom_types:
                    # Atom type already exists in morphology, so increment the
                    # atom_type number by one
                    for i in range(1, len(atom_type)):
                        # Work out where the integer start so we can increment
                        # it (should be i = 1 for one-character element names)
                        try:
                            integer = int(atom_type[i:])
                            break
                        except:
                            continue
                    atom_type = "".join([atom_type[:i], str(integer + 1)])
                morphology_atom_types.append(atom_type)
                mapping_for_this_FF[lj_interaction[0]] = atom_type
            force_field_mappings.append(mapping_for_this_FF)
        for CG_site, directory in CG_to_template_dirs.items():
            FF_loc = os.path.join(directory, CG_to_template_force_fields[CG_site])
            CG_to_template_mappings[CG_site] = force_field_mappings[
                force_field_locations.index(FF_loc)
            ]
        return CG_to_template_mappings


class atomistic:
    def __init__(
        self,
        molecule_index,
        site_IDs,
        CG_dictionary,
        molecule_lengths,
        rolling_AA_index,
        ghost_dictionary,
        parameter_dict,
    ):
        # This class sees individual molecules.
        self.no_atoms_in_morphology = rolling_AA_index
        self.molecule_index = molecule_index
        self.molecule_lengths = molecule_lengths
        self.site_IDs = site_IDs
        self.CG_dictionary = CG_dictionary
        # Get the dictionary of all the CG sites in this molecule
        # self.CGMonomerDictionary = self.getCGMonomerDict()
        # Import the parXX.py parameters
        for key, value in parameter_dict.items():
            self.__dict__[key] = value
        self.AA_templates_dictionary = {}
        # Load the template file for each CG atom
        for CG_atom_type in list(self.CG_to_template_files.keys()):
            template_dictionary = hf.load_morphology_xml(
                os.path.join(
                    self.CG_to_template_dirs[CG_atom_type],
                    self.CG_to_template_files[CG_atom_type],
                )
            )
            template_dictionary = self.remap_atom_types(
                template_dictionary, parameter_dict["new_type_mappings"][CG_atom_type]
            )
            template_dictionary = hf.add_unwrapped_positions(template_dictionary)
            self.AA_templates_dictionary[CG_atom_type] = template_dictionary
        fine_grained = self.run_fine_grainer(ghost_dictionary)
        self.AA_dictionary = fine_grained[0]
        self.atom_ID_lookup_table = fine_grained[1]
        self.ghost_dictionary = fine_grained[2]

    def return_data(self):
        # Return the important fine-grained results from this class
        return self.AA_dictionary, self.atom_ID_lookup_table, self.ghost_dictionary

    def get_CG_monomer_dict(self):
        CG_monomer_dictionary = {
            "position": [],
            "image": [],
            "mass": [],
            "diameter": [],
            "type": [],
            "body": [],
            "bond": [],
            "angle": [],
            "dihedral": [],
            "improper": [],
            "charge": [],
            "lx": 0,
            "ly": 0,
            "lz": 0,
            "xy": 0,
            "xz": 0,
            "yz": 0,
        }
        # First, do just the positions and find the newsiteIDs for each CG site
        for site_ID in self.site_IDs:
            CG_monomer_dictionary["position"].append(
                self.CG_dictionary["position"][site_ID]
            )
        # Now sort out the other one-per-atom properties
        for key in ["image", "mass", "diameter", "type", "body", "charge"]:
            if len(self.CG_dictionary[key]) != 0:
                for site_ID in self.site_IDs:
                    CG_monomer_dictionary[key].append(self.CG_dictionary[key][site_ID])
        # Now rewrite the bonds based on the newsiteIDs
        for key in ["bond", "angle", "dihedral", "improper"]:
            for element in self.CG_dictionary[key]:
                for site_ID in self.site_IDs:
                    if (site_ID in element) and (
                        element not in CG_monomer_dictionary[key]
                    ):
                        CG_monomer_dictionary[key].append(element)
        # Now update the box parameters
        for key in ["lx", "ly", "lz"]:
            CG_monomer_dictionary[key] = self.CG_dictionary[key]
        CG_monomer_dictionary = hf.add_unwrapped_positions(CG_monomer_dictionary)
        CG_monomer_dictionary["natoms"] = len(CG_monomer_dictionary["position"])
        return CG_monomer_dictionary

    def remap_atom_types(self, template_dict, mappings):
        # Remap the atom types first
        for index, atom_type in enumerate(template_dict["type"]):
            template_dict["type"][index] = mappings[atom_type]
        # Then rename the constraints appropriately
        constraint_types = ["bond", "angle", "dihedral", "improper"]
        for constraint_type in constraint_types:
            for constraint_index, constraint in enumerate(
                template_dict[constraint_type]
            ):
                new_constraint0 = "-".join(
                    [template_dict["type"][atom_ID] for atom_ID in constraint[1:]]
                )
                # Assign the constraint label
                template_dict[constraint_type][constraint_index][0] = new_constraint0
        return template_dict

    def run_fine_grainer(self, ghost_dictionary):
        AA_dictionary = {
            "position": [],
            "image": [],
            "unwrapped_position": [],
            "mass": [],
            "diameter": [],
            "type": [],
            "body": [],
            "bond": [],
            "angle": [],
            "dihedral": [],
            "improper": [],
            "charge": [],
            "lx": 0,
            "ly": 0,
            "lz": 0,
            "xy": 0,
            "xz": 0,
            "yz": 0,
        }
        # Find the COMs of each CG site in the system, so that we know where to
        # move the template to
        CG_co_ms, self.CG_to_template_AAIDs = self.get_AA_template_position(
            self.CG_to_template_AAIDs
        )
        # Need to keep track of the atom ID numbers globally - run_fine_grainer
        # sees individual monomers, atomistic sees molecules and the xml needs
        # to contain the entire morphology.
        no_atoms_in_molecule = 0
        CG_type_list = {}
        for site_ID in self.site_IDs:
            CG_type_list[site_ID] = self.CG_dictionary["type"][site_ID]
        # Sort the CG sites into monomers so we can iterate over each monomer
        # in order to perform the fine-graining
        monomer_list = self.sort_into_monomers(CG_type_list)
        current_monomer_index = sum(self.molecule_lengths[: self.molecule_index])
        atom_ID_lookup_table = {}
        # Calculate the total number of permitted atoms
        total_permitted_atoms = self.get_total_permitted_atoms(monomer_list)
        # Set the initial and final atom indices to None initially, so that we
        # don't add terminating units for small molecules
        start_atom_index = None
        end_atom_index = None
        for monomer_no, monomer in enumerate(monomer_list):
            # This monomer should have the same template file for all CG sites
            # in the monomer, if not we've made a mistake in splitting the
            # monomers. So check this:
            template_files = []
            monomer_CG_types = []
            for CG_site in monomer:
                template_files.append(
                    self.CG_to_template_files[self.CG_dictionary["type"][CG_site]]
                )
                monomer_CG_types.append(self.CG_dictionary["type"][CG_site])
            if len(set(template_files)) != 1:
                print(monomer)
                print(monomerCGTypes)
                print(templateFiles)
                raise SystemError("Not all monomer sites are the same template")
            # Copy the template dictionary for editing for this monomer
            this_monomer_dictionary = copy.deepcopy(
                self.AA_templates_dictionary[self.CG_dictionary["type"][monomer[0]]]
            )
            for key in ["lx", "ly", "lz"]:  # TODO: Tilts
                this_monomer_dictionary[key] = self.CG_dictionary[key]
            # Include the image tag in case it's not present in the template
            if len(this_monomer_dictionary["image"]) == 0:
                this_monomer_dictionary["image"] = [[0, 0, 0]] * len(
                    this_monomer_dictionary["position"]
                )
            for site_ID in monomer:
                site_posn = np.array(self.CG_dictionary["unwrapped_position"][site_ID])
                site_translation = site_posn - CG_co_ms[CG_type_list[site_ID]]
                atom_ID_lookup_table[site_ID] = [
                    CG_type_list[site_ID],
                    [
                        x + no_atoms_in_molecule + self.no_atoms_in_morphology
                        for x in self.CG_to_template_AAIDs[CG_type_list[site_ID]]
                    ],
                ]
                # Add the atoms in based on the CG site position
                for AAID in self.CG_to_template_AAIDs[CG_type_list[site_ID]]:
                    this_monomer_dictionary["unwrapped_position"][AAID] = list(
                        np.array(this_monomer_dictionary["unwrapped_position"][AAID])
                        + site_translation
                    )
                # Next sort out the rigid bodies
                if CG_type_list[site_ID] in self.rigid_body_sites:
                    # Every rigid body needs a ghost particle that describes
                    # its CoM
                    AAID_positions = []
                    AAID_atom_types = []
                    # If the key is specified with no values, assume that all
                    # the AAIDs in the template constitute the rigid body
                    if len(self.rigid_body_sites[CG_type_list[site_ID]]) == 0:
                        self.rigid_body_sites[CG_type_list[site_ID]] = list(
                            np.arange(
                                len(self.CG_to_template_AAIDs[CG_type_list[site_ID]])
                            )
                        )
                    for AAID in self.rigid_body_sites[CG_type_list[site_ID]]:
                        this_monomer_dictionary["body"][AAID] = current_monomer_index
                        AAID_positions.append(
                            this_monomer_dictionary["unwrapped_position"][AAID]
                        )
                        AAID_atom_types.append(this_monomer_dictionary["type"][AAID])
                    # Now create the ghost particle describing the rigid body
                    ghost_COM = hf.calc_COM(
                        AAID_positions, list_of_atom_types=AAID_atom_types
                    )
                    ghost_dictionary["unwrapped_position"].append(ghost_COM)
                    ghost_dictionary["mass"].append(1.0)
                    ghost_dictionary["diameter"].append(1.0)
                    ghost_dictionary["type"].append(
                        "R{:s}".format(CG_type_list[site_ID])
                    )
                    ghost_dictionary["body"].append(current_monomer_index)
                    ghost_dictionary["charge"].append(0.0)
                    # Then create the corresponding CG anchorpoint
                    ghost_dictionary["unwrapped_position"].append(ghost_COM)
                    ghost_dictionary["mass"].append(1.0)
                    ghost_dictionary["diameter"].append(1.0)
                    ghost_dictionary["type"].append(
                        "X{:s}".format(CG_type_list[site_ID])
                    )
                    ghost_dictionary["body"].append(-1)
                    ghost_dictionary["charge"].append(0.0)
                    # Now create a bond between them
                    # We want to bond together the previous two ghost particles,
                    # so this should work as it requires no knowledge of the
                    # number of ghost particles already in the system.
                    ghost_dictionary["bond"].append(
                        [
                            "{0:s}-{1:s}".format(
                                ghost_dictionary["type"][-2],
                                ghost_dictionary["type"][-1],
                            ),
                            "".join(["*" + str(len(ghost_dictionary["type"]) - 2)]),
                            "".join(["*" + str(len(ghost_dictionary["type"]) - 1)]),
                        ]
                    )
                else:
                    # Create a ghost particle that describe the CG anchorpoint
                    # for the non-rigid body
                    ghost_dictionary["unwrapped_position"].append(
                        self.CG_dictionary["unwrapped_position"][site_ID]
                    )
                    ghost_dictionary["mass"].append(1.0)
                    ghost_dictionary["diameter"].append(1.0)
                    ghost_dictionary["type"].append(
                        "X{:s}".format(CG_type_list[site_ID])
                    )
                    ghost_dictionary["body"].append(-1)
                    ghost_dictionary["charge"].append(0.0)
                    # Add in bonds between the CG anchorpoints and the atom
                    # closest to the CG site.
                    # Find the atom closest to the CG site
                    closest_atom_ID = None
                    closest_atom_posn = 1E99
                    for AAID, AA_position in enumerate(
                        this_monomer_dictionary["unwrapped_position"]
                    ):
                        separation = hf.calculate_separation(
                            self.CG_dictionary["unwrapped_position"][site_ID],
                            AA_position,
                        )
                        if separation < closest_atom_posn:
                            closest_atom_posn = separation
                            closest_atom_ID = AAID
                    # Add in the bond:
                    # Note that, in order to distinguish between the
                    # ghost_atom_IDs and the real_atomIDs, I've put an
                    # underscore in front of the closest_atom_ID, and a * in
                    # front of the ghost_atom_ID. When incrementing the
                    # atom_IDs for this monomer or this molecule, the
                    # real_atom_IDs will be incremented correctly.
                    # Later, when the ghost atoms are added to the main system,
                    # the ghost_atom_IDs will be incremented according to the
                    # number of atoms in the whole system (i.e. the ghosts
                    # appear at the end of the real atoms).
                    # At this time, the real_atom_IDs will be left unchanged
                    # because they are already correct for the whole system.
                    ghost_dictionary["bond"].append(
                        [
                            "{0:s}-{1:s}".format(
                                ghost_dictionary["type"][-1],
                                this_monomer_dictionary["type"][closest_atom_ID],
                            ),
                            "".join(["*" + str(len(ghost_dictionary["type"]) - 1)]),
                            "".join(
                                ["_" + str(closest_atom_ID + no_atoms_in_molecule)]
                            ),
                        ]
                    )
            # Now add in the bonds between CG_sites in this monomer
            for bond in self.CG_dictionary["bond"]:
                if (bond[1] in monomer) and (bond[2] in monomer):
                    CG_bond_type = bond[0]
                    this_monomer_dictionary["bond"].append(
                        self.CG_to_template_bonds[CG_bond_type]
                    )
            # Now need to add in the additional_constraints for this monomer
            # (which include the bond, angle and dihedral for the inter-monomer
            # connections). However, we need a check to make sure that we don't
            # add stuff for the final monomer (because those atoms +25 don't
            # exist in this molecule!)
            for constraint in self.additional_constraints:
                # Firstly, skip this constraint if the current monomer doesn't
                # have the relevant atom types
                if (
                    all(
                        [
                            atom_type in set(this_monomer_dictionary["type"])
                            for atom_type in constraint[0].split("-")
                        ]
                    )
                    is False
                ):
                    continue
                # Also check that we're not at the final monomer
                at_final_monomer = False
                for atom_ID in constraint[2:]:
                    if (no_atoms_in_molecule + atom_ID + 1) > total_permitted_atoms:
                        at_final_monomer = True
                        break
                if at_final_monomer is True:
                    break
                # Work out which key to write the constraint to based on its
                # length:
                # 3 = Bond, 4 = Angle, 5 = Dihedral, 6 = Improper
                constraint_type = ["bond", "angle", "dihedral", "improper"]
                this_monomer_dictionary[constraint_type[len(constraint) - 3]].append(
                    constraint
                )
            # Finally, increment the atom IDs to take into account previous
            # monomers in this molecule and then update the AADictionary.
            # Note that the ghost dictionary bond was already updated to have
            # the correct realAtom AAID for this molecule when the bond was
            # created. Therefore, leave the ghost dictionary unchanged.
            this_monomer_dictionary, ghost_dictionary = hf.increment_atom_IDs(
                this_monomer_dictionary,
                ghost_dictionary,
                no_atoms_in_molecule,
                modify_ghost_dictionary=False,
            )
            no_atoms_in_molecule += len(this_monomer_dictionary["type"])
            current_monomer_index += 1
            # Update the current AA dictionary with this monomer
            AA_dictionary = self.update_molecule_dictionary(
                this_monomer_dictionary, AA_dictionary
            )
        # All Monomers sorted, now for the final bits
        AA_dictionary["natoms"] = no_atoms_in_molecule
        for key in ["lx", "ly", "lz"]:
            AA_dictionary[key] = this_monomer_dictionary[key]
        # Now we terminate the molecules using the technique in the
        # add_hydrogens_to_UA analysis script
        new_hydrogen_data = []
        for atom_index, atom_type in enumerate(AA_dictionary["type"]):
            if atom_type not in self.molecule_terminating_connections.keys():
                continue
            bonded_AAIDs = []
            # Iterate over all termination connections defined for this
            # atom_type (in case we are trying to do something mega
            # complicated)
            for connection_info in self.molecule_terminating_connections[atom_type]:
                for [bond_name, AAID1, AAID2] in AA_dictionary["bond"]:
                    if AAID1 == atom_index:
                        if AAID2 not in bonded_AAIDs:
                            bonded_AAIDs.append(AAID2)
                    elif AAID2 == atom_index:
                        if AAID1 not in bonded_AAIDs:
                            bonded_AAIDs.append(AAID1)
                if len(bonded_AAIDs) != connection_info[0]:
                    continue
                new_hydrogen_positions = hf.get_terminating_positions(
                    AA_dictionary["unwrapped_position"][atom_index],
                    [
                        AA_dictionary["unwrapped_position"][bonded_AAID]
                        for bonded_AAID in bonded_AAIDs
                    ],
                    1,
                )
                for hydrogen_position in new_hydrogen_positions:
                    new_hydrogen_data.append([atom_index, list(hydrogen_position)])
        AA_dictionary = self.add_terminating_to_molecule(
            AA_dictionary, new_hydrogen_data
        )
        AA_dictionary = hf.add_wrapped_positions(AA_dictionary)
        AA_dictionary = hf.add_masses(AA_dictionary)
        AA_dictionary = hf.add_diameters(AA_dictionary)
        # Now the molecule is done, we need to add on the correct identifying
        # numbers for all the bonds, angles and dihedrals (just as we did
        # between monomers) for the other molecules in the system, so that they
        # all connect to the right atoms.
        # Note that here we need to increment the '_'+ATOMIDs in the ghost
        # dictionary to take into account the number of molecules.
        AA_dictionary, ghost_dictionary = hf.increment_atom_IDs(
            AA_dictionary,
            ghost_dictionary,
            self.no_atoms_in_morphology,
            modify_ghost_dictionary=True,
        )
        return AA_dictionary, atom_ID_lookup_table, ghost_dictionary

    def add_terminating_to_molecule(self, AA_dictionary, new_hydrogen_data):
        hydrogen_ID = len(AA_dictionary["type"])
        for hydrogen in new_hydrogen_data:
            # hydrogen has the form [BondedAtomID, [position]]
            AA_dictionary["unwrapped_position"].append([coord for coord in hydrogen[1]])
            AA_dictionary["type"].append("H1")
            AA_dictionary["body"].append(-1)
            AA_dictionary["charge"].append(0.0)
            AA_dictionary["bond"].append(
                [
                    "-".join([AA_dictionary["type"][hydrogen[0]], "h1"]),
                    hydrogen[0],
                    hydrogen_ID,
                ]
            )
            AA_dictionary["natoms"] += 1
            hydrogen_ID += 1
        return AA_dictionary

    def get_total_permitted_atoms(self, monomer_list):
        # Work out how many atoms we have in the molecule so that we don't
        # create any constraints including atoms outside this molecule
        total_permitted_atoms = 0
        for monomer in monomer_list:
            for CG_site_ID in monomer:
                CG_site_type = self.CG_dictionary["type"][CG_site_ID]
                total_permitted_atoms += len(self.CG_to_template_AAIDs[CG_site_type])
        return total_permitted_atoms

    def sort_into_monomers(self, type_list_sequence):
        monomer_list = []
        molecule_site_IDs = copy.deepcopy(self.site_IDs)
        # Iterate over the entire bondlist until it's done
        bond_list = copy.deepcopy(self.CG_dictionary["bond"])
        while len(molecule_site_IDs) > 0:
            # Add the first atom to the first monomer
            this_monomer = []
            monomer_types_added = [type_list_sequence[molecule_site_IDs[0]]]
            this_monomer.append(molecule_site_IDs[0])
            added_new_site = True
            while added_new_site is True:
                # Keep adding new, connected atoms until we can't add any more
                added_new_site = False
                bond_pop_list = []
                # Find bonded atoms that are not of the same type
                for bond_no, bond in enumerate(bond_list):
                    if (
                        (bond[1] in this_monomer)
                        and (bond[2] not in this_monomer)
                        and (type_list_sequence[bond[2]] not in monomer_types_added)
                    ):
                        this_monomer.append(bond[2])
                        monomer_types_added.append(type_list_sequence[bond[2]])
                        added_new_site = True
                    elif (
                        (bond[2] in this_monomer)
                        and (bond[1] not in this_monomer)
                        and (type_list_sequence[bond[1]] not in monomer_types_added)
                    ):
                        this_monomer.append(bond[1])
                        monomer_types_added.append(type_list_sequence[bond[1]])
                        added_new_site = True
                    elif (bond[1] in this_monomer) and (bond[2] in this_monomer):
                        pass
                    else:
                        continue
                    bond_pop_list.append(bond_no)
                # Remove the bonds we've already considered
                bond_pop_list.sort(reverse=True)
                for bond_index in bond_pop_list:
                    bond_list.pop(bond_index)
            for atom_index in this_monomer:
                molecule_site_IDs.remove(atom_index)
            monomer_list.append(this_monomer)
        monomer_types_list = []
        for monomer in monomer_list:
            monomer_types_list.append([])
            for atom in monomer:
                monomer_types_list[-1].append(type_list_sequence[atom])
        return monomer_list

    def update_molecule_dictionary(self, current_monomer_dictionary, AA_dictionary):
        # Update AA_dictionary with all of the values in
        # current_monomer_dictionary, except ths system dimensions which will
        # be sorted later
        key_list = list(AA_dictionary.keys())
        key_list.remove("lx")
        key_list.remove("ly")
        key_list.remove("lz")
        key_list.remove("xy")
        key_list.remove("xz")
        key_list.remove("yz")
        for key in key_list:
            for value in current_monomer_dictionary[key]:
                AA_dictionary[key].append(value)
        return AA_dictionary

    def get_AA_template_position(self, CG_to_template_AAIDs):
        CG_COMs = {}
        # For each CG site, determine the types and positions so we can
        # calculate the COM
        for site_name in list(CG_to_template_AAIDs.keys()):
            atom_IDs = CG_to_template_AAIDs[site_name]
            AA_template = self.AA_templates_dictionary[site_name]
            # If the key's length is zero, then add all the atoms from the
            # template
            if len(atom_IDs) == 0:
                atom_IDs = list(np.arange(len(AA_template["type"])))
                # Update self.CGToTemplateAAIDs with these for later on
                CG_to_template_AAIDs[site_name] = atom_IDs
            site_positions = []
            site_types = []
            for atom_ID in atom_IDs:
                site_types.append(AA_template["type"][atom_ID])
                site_positions.append(AA_template["unwrapped_position"][atom_ID])
            # These output as numpy arrays because we can't do maths with lists
            CG_COMs[site_name] = hf.calc_COM(
                site_positions, list_of_atom_types=site_types
            )
        return CG_COMs, CG_to_template_AAIDs
