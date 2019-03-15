import os
import sys
import numpy as np
from morphct.code import helper_functions as hf


def check_bonds(morphology, bond_dict):
    for bond in morphology["bond"]:
        posn1 = np.array(morphology["position"][bond[1]])
        +(
            np.array(morphology["image"][bond[1]])
            * np.array([morphology["lx"], morphology["ly"], morphology["lz"]])
        )
        posn2 = np.array(morphology["position"][bond[2]])
        +(
            np.array(morphology["image"][bond[2]])
            * np.array([morphology["lx"], morphology["ly"], morphology["lz"]])
        )
        separation = hf.calculate_separation(posn1, posn2)
        if separation >= morphology["lx"] / 2.0:
            print(
                "Periodic bond found:",
                bond,
                "because separation =",
                separation,
                ">=",
                morphology["lx"] / 2.0,
            )
            morphology = move_bonded_atoms(bond[1], morphology, bond_dict)
    return morphology


def zero_out_images(morphology):
    morphology["image"] = [[0, 0, 0]] * len(morphology["position"])
    return morphology


def get_bond_dict(morphology):
    bond_dict = {atom_id: [] for atom_id, atom_type in enumerate(morphology["type"])}
    for bond in morphology["bond"]:
        bond_dict[bond[1]].append(bond[2])
        bond_dict[bond[2]].append(bond[1])
    return bond_dict


def move_bonded_atoms(central_atom, morphology, bond_dict):
    for bonded_atom in bond_dict[central_atom]:
        atom1_posn = morphology["position"][central_atom]
        atom2_posn = morphology["position"][bonded_atom]
        sep_vec = np.array(atom1_posn) - np.array(atom2_posn)
        moved = False
        for axis, value in enumerate(sep_vec):
            if value > morphology["lx"] / 2.0:
                morphology["position"][bonded_atom][axis] += morphology["lx"]
                moved = True
            if value < -morphology["lx"] / 2.0:
                morphology["position"][bonded_atom][axis] -= morphology["lx"]
                moved = True
        if moved:
            morphology = move_bonded_atoms(bonded_atom, morphology, bond_dict)
    return morphology


def main():
    list_of_files = sys.argv[1:]
    if len(list_of_files) < 1:
        print("No files requested to convert!")
        exit()
    for file_name in list_of_files:
        print("Fixing the images for {:s}...".format(file_name))
        morphology = hf.load_morphology_xml(file_name)
        morphology = zero_out_images(morphology)
        bond_dict = get_bond_dict(morphology)
        morphology = check_bonds(morphology, bond_dict)
        file_directory, split_file_name = os.path.split(file_name)
        hf.write_morphology_xml(
            morphology,
            os.path.join(file_directory, "".join(["image_fix_", split_file_name])),
        )


if __name__ == "__main__":
    main()
