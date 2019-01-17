import numpy as np
from morphct.code import helper_functions as hf
import argparse


def calculate_hydrogen_positions(morphology_dict, hydrogens_to_add):
    """This function calculates the position of the hydrogen based
    on the number and positions of the other bonded species, and
    the number of hydrogens required to be added to the atom"""
    hydrogen_positions = []
    # First create a lookup table that describes exactly how many bonds
    # each atom has to neighbours
    number_of_bonds = {x: [0, []] for x in range(len(morphology_dict["type"]))}
    # Keys == atom_ID, values == [number_of_bonded_species, [bonded_atom_IDs]]
    for bond in morphology_dict["bond"]:
        atom_ID1 = bond[1]
        atom_ID2 = bond[2]
        if atom_ID2 not in number_of_bonds[atom_ID1][1]:
            number_of_bonds[atom_ID1][0] += 1
            number_of_bonds[atom_ID1][1].append(atom_ID2)
        if atom_ID1 not in number_of_bonds[atom_ID2][1]:
            number_of_bonds[atom_ID2][0] += 1
            number_of_bonds[atom_ID2][1].append(atom_ID1)
    for atom_ID, atom_type in enumerate(morphology_dict["type"]):
        # Skip if we don't have to add hydrogens to the current atom's type
        if atom_type not in hydrogens_to_add.keys():
            continue
        for bond_definition in hydrogens_to_add[atom_type]:
            # Skip if the current atom does not have the right number of bonds
            if number_of_bonds[atom_ID][0] != bond_definition[0]:
                continue
            # Otherwise, we need to add hydrogens_to_add[atom_type][1] hydrogens
            # to this atom
            current_atom_posn = morphology_dict["unwrapped_position"][atom_ID]
            # First get the vector to the average position of the bonded
            # neighbours
            average_position_of_bonded_atoms = np.array([0.0, 0.0, 0.0])
            for bonded_atom in number_of_bonds[atom_ID][1]:
                bond_vector = (
                    np.array(morphology_dict["unwrapped_position"][bonded_atom])
                    - current_atom_posn
                )
                bond_vector /= np.linalg.norm(bond_vector)
                average_position_of_bonded_atoms += bond_vector
            [x, y, z] = current_atom_posn + (
                -1.06
                * (
                    average_position_of_bonded_atoms
                    / np.linalg.norm(average_position_of_bonded_atoms)
                )
            )
            if bond_definition[1] == 1:
                # Easy, this is the perylene code
                # Simply reverse the bonded vector and make it the hydrogen
                # position at a distance of 1.06 angstroems
                hydrogen_positions.append([int(atom_ID), np.array([x, y, z])])
            # Initial position for all hydrogens
            elif bond_definition[1] == 2:
                # As above (to get the right plane), but then rotated
                # +(109.5/2) degrees and -(109.5/2) degrees around the bonding
                # axis
                rotation_axis = np.array(
                    morphology_dict["unwrapped_position"][
                        number_of_bonds[atom_ID][1][0]
                    ]
                ) - np.array(
                    morphology_dict["unwrapped_position"][
                        number_of_bonds[atom_ID][1][-1]
                    ]
                )
                rotation_axis /= np.linalg.norm(rotation_axis)
                # Rotation matrix calculations from:
                # http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
                # The array that describes the 3D rotation of (x, y, z) around
                # the point (a, b, c) through
                # the unit axis <u, v, w> by the angle theta is given by:
                # [ (a(v^2 + w^2) - u(bv + cw - ux - vy - wz))(1 - cos(theta))
                #       + x*cos(theta) + (-cv + bw - wy + vz)sin(theta),
                #   (b(u^2 + w^2) - v(au + cw - ux - vy - wz))(1 - cos(theta))
                #       + y*cos(theta) + (cu - aw + wx - uz)sin(theta),
                #   (c(u^2 + v^2) - w(au + bv - ux - vy - wz))(1 - cos(theta))
                #       + z*cos(theta) + (-bu + av - vx + uy)sin(theta) ]
                [a, b, c] = current_atom_posn
                [u, v, w] = rotation_axis
                for theta in [
                    (109.5 / 2.0) * (np.pi / 180.0),
                    -(109.5 / 2.0) * (np.pi / 180.0),
                ]:
                    new_position = np.array(
                        [
                            (
                                a * (v ** 2 + w ** 2)
                                - u * ((b * v) + (c * w) - (u * x) - (v * y) - (w * z))
                            )
                            * (1 - np.cos(theta))
                            + (x * np.cos(theta))
                            + (
                                (-(c * v) + (b * w) - (w * y) + (v * z)) * np.sin(theta)
                            ),
                            (
                                b * (u ** 2 + w ** 2)
                                - v * ((a * u) + (c * w) - (u * x) - (v * y) - (w * z))
                            )
                            * (1 - np.cos(theta))
                            + (y * np.cos(theta))
                            + (((c * u) - (a * w) + (w * x) - (u * z)) * np.sin(theta)),
                            (
                                c * (u ** 2 + v ** 2)
                                - w * ((a * u) + (b * v) - (u * x) - (v * y) - (w * z))
                            )
                            * (1 - np.cos(theta))
                            + (z * np.cos(theta))
                            + (
                                (-(b * u) + (a * v) - (v * x) + (u * y)) * np.sin(theta)
                            ),
                        ]
                    )
                    hydrogen_positions.append([int(atom_ID), new_position])
            elif bond_definition[1] == 3:
                # As for one (to get the right side of the bonded atom),
                # rotate the first one up by 70.5 (180 - 109.5) and then rotate
                # around by 109.5 degrees for the other two.
                # The first hydrogen can be rotated around any axis
                # perpendicular to the only bond present
                axis_to_bond = current_atom_posn - np.array(
                    morphology_dict["unwrapped_position"][
                        number_of_bonds[atom_ID][1][0]
                    ]
                )
                # Now find one of the set of vectors [i, j, k] perpendicular to
                # this one so we can place the first hydrogen.
                # Do this by setting i = j = 1 and solve for k (given that
                # current_atom_posn[0]*i + current_atom_posn[1]*j
                # + current_atom_posn[2]*k = 0)
                first_hydrogen_rotation_axis = np.array(
                    [1, 1, -(axis_to_bond[0] + axis_to_bond[1]) / axis_to_bond[2]]
                )
                first_hydrogen_rotation_axis /= np.linalg.norm(
                    first_hydrogen_rotation_axis
                )

                [a, b, c] = current_atom_posn
                [u, v, w] = first_hydrogen_rotation_axis
                # First hydrogen
                theta = 70.5 * np.pi / 180.0
                new_position = np.array(
                    [
                        (
                            a * (v ** 2 + w ** 2)
                            - u * ((b * v) + (c * w) - (u * x) - (v * y) - (w * z))
                        )
                        * (1 - np.cos(theta))
                        + (x * np.cos(theta))
                        + ((-(c * v) + (b * w) - (w * y) + (v * z)) * np.sin(theta)),
                        (
                            b * (u ** 2 + w ** 2)
                            - v * ((a * u) + (c * w) - (u * x) - (v * y) - (w * z))
                        )
                        * (1 - np.cos(theta))
                        + (y * np.cos(theta))
                        + (((c * u) - (a * w) + (w * x) - (u * z)) * np.sin(theta)),
                        (
                            c * (u ** 2 + v ** 2)
                            - w * ((a * u) + (b * v) - (u * x) - (v * y) - (w * z))
                        )
                        * (1 - np.cos(theta))
                        + (z * np.cos(theta))
                        + ((-(b * u) + (a * v) - (v * x) + (u * y)) * np.sin(theta)),
                    ]
                )
                hydrogen_positions.append([int(atom_ID), new_position])
                # Second and third hydrogens
                # Rotate these from the newPosition +/-120 degrees around the
                # vector axisToBond from the position currentAtomPosn -
                # axisToBond
                [x, y, z] = new_position
                [a, b, c] = current_atom_posn + (np.cos(theta) * axis_to_bond)
                [u, v, w] = (np.cos(theta) * axis_to_bond) / np.linalg.norm(
                    np.cos(theta) * axis_to_bond
                )
                for theta in [120 * (np.pi / 180.0), -120 * (np.pi / 180.0)]:
                    new_hydrogen_position = np.array(
                        [
                            (
                                a * (v ** 2 + w ** 2)
                                - u * ((b * v) + (c * w) - (u * x) - (v * y) - (w * z))
                            )
                            * (1 - np.cos(theta))
                            + (x * np.cos(theta))
                            + (
                                (-(c * v) + (b * w) - (w * y) + (v * z)) * np.sin(theta)
                            ),
                            (
                                b * (u ** 2 + w ** 2)
                                - v * ((a * u) + (c * w) - (u * x) - (v * y) - (w * z))
                            )
                            * (1 - np.cos(theta))
                            + (y * np.cos(theta))
                            + (((c * u) - (a * w) + (w * x) - (u * z)) * np.sin(theta)),
                            (
                                c * (u ** 2 + v ** 2)
                                - w * ((a * u) + (b * v) - (u * x) - (v * y) - (w * z))
                            )
                            * (1 - np.cos(theta))
                            + (z * np.cos(theta))
                            + (
                                (-(b * u) + (a * v) - (v * x) + (u * y)) * np.sin(theta)
                            ),
                        ]
                    )
                    hydrogen_positions.append([int(atom_ID), new_hydrogen_position])
    return hydrogen_positions


def add_hydrogens_to_morph(morphology_dict, hydrogen_positions):
    """This function adds the hydrogen atoms into the morphology
    to be exported as the AA xml"""
    # hydrogen_positions are of the following format:
    # [carbon_index, position_of_hydrogen]
    for hydrogen_atom in hydrogen_positions:
        morphology_dict["type"].append("H1")
        morphology_dict["unwrapped_position"].append(hydrogen_atom[1])
        morphology_dict["bond"].append(
            ["C-H", hydrogen_atom[0], morphology_dict["natoms"]]
        )
        morphology_dict["natoms"] += 1
        morphology_dict["mass"].append(1.00794)
        morphology_dict["charge"].append(0.0)
        morphology_dict["diameter"].append(0.53)
        morphology_dict["body"].append(morphology_dict["body"][hydrogen_atom[0]])
    return morphology_dict


def find_information(args):
    """
    Function determines the dictionary and sigma
    values for the given argument value (or asks
    for one.)
    Requires:
        None
    Returns:
        Dictionary of atom types and number of hydrogens
        sigma
    """

    built_ins = {
        "PCBM": [{"CHA": [[2, 1]], "CH2": [[2, 2]], "CE": [[1, 3]]}, 1.0],
        "P3HT": [{"CA": [[2, 1]], "CT": [[2, 2], [1, 3]]}, 3.905],
        "P3HT_sig_1": [{"CA": [[2, 1]], "CT": [[2, 2], [1, 3]]}, 1.0],
        "DBP": [{"C": [[2, 1]], "CA": [[2, 1]], "CT": [[2, 2], [1, 3]]}, 3.905],
        "PERYLENE": [{"C": [[2, 1]]}, 3.8],
        "BDT-TPD": [
            {
                "CS": [[2, 1]],
                "C!": [[2, 1]],
                "C*": [[2, 1]],
                "CT": [[2, 2], [1, 3], [3, 1]],
                "CP": [[2, 1]],
            },
            3.55,
        ],
        "BDT-TPD/PCBM": [
            {
                "CS": [[2, 1]],
                "C!": [[2, 1]],
                "CT": [[2, 2], [1, 3], [3, 1]],
                "CP": [[2, 1]],
                "FCT": [[2, 2], [1, 3]],
            },
            3.55,
        ],
        "P3HT/PCBM": [
            {
                "CA": [[2, 1]],
                "CT": [[2, 2], [1, 3]],
                "FCA": [[2, 1]],
                "FCT": [[2, 2], [1, 3]],
            },
            3.905,
        ],
    }

    if args.molecule_source is None:
        print("You have not specified the dictionary you want to use to add hydrogens.")
        print("Current built-ins are:")
        for key, val in built_ins.items():
            print(key)
        print("Or you can specify a python script with the dictionary.")
        system_name = input("Which would you like? >> ")
        if system_name == "":
            exit("No information given. Exiting.")
    else:
        system_name = args.molecule_source
    if len(system_name.split(".")) > 1:
        # Import python script
        file_to_load = system_name.split(".")[0]
        new_module = __import__(file_to_load)
        return new_module.DataToLoad()
    else:
        return built_ins[system_name][0], built_ins[system_name][1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m",
        "--molecule_source",
        required=False,
        help="""Specify the source of the add_hydrogens dictionary and sigma value.
        This can be in the form of a python script which returns a function of the
        dictionary and sigma or the name of the molecule.""",
    )
    args, input_files = parser.parse_known_args()
    hydrogens_to_add, sigma_val = find_information(args)
    for input_file in input_files:
        # This dictionary has keys of the atom type, and values where the first
        # element is the number of bonds required for us to add a hydrogen to
        # the atom and the second element of the value defines how many
        # hydrogens to add to said atom.
        print(
            "THIS FUNCTION IS SET UP TO USE A DICTIONARY TO DEFINE HOW MANY HYDROGENS "
            " TO ADD TO BONDS OF A SPECIFIC TYPE WITH A CERTAIN NUMBER OF BONDS"
        )
        print(hydrogens_to_add)
        print(
            "IF THE ABOVE DICTIONARY DOESN'T LOOK RIGHT, PLEASE TERMINATE NOW AND "
            " IGNORE ANY OUTPUTS UNTIL THE DICTIONARY HAS BEEN RECTIFIED"
        )
        print("Additionally, we're using a sigma value of", sigma_val)
        morphology_dict = hf.load_morphology_xml(input_file, sigma=sigma_val)
        morphology_dict = hf.add_unwrapped_positions(morphology_dict)
        hydrogen_positions = calculate_hydrogen_positions(
            morphology_dict, hydrogens_to_add
        )
        morphology_dict = add_hydrogens_to_morph(morphology_dict, hydrogen_positions)
        morphology_dict = hf.add_wrapped_positions(morphology_dict)
        hf.write_morphology_xml(
            morphology_dict,
            input_file.replace(".xml", "_AA.xml"),
            check_wrapped_posns=False,
        )


if __name__ == "__main__":
    main()
