import copy
import csv
import os
import pickle
import sys
import multiprocessing as mp
import numpy as np
from morphct.definitions import PROJECT_ROOT


try:
    from hoomd_script import init
    has_hoomd = True
except ImportError:
    has_hoomd = False


# UNIVERSAL CONSTANTS, DO NOT CHANGE!
elementary_charge = 1.60217657E-19  # C
k_B = 1.3806488E-23  # m^{2} kg s^{-2} K^{-1}
hbar = 1.05457173E-34  # m^{2} kg s^{-1}


sys.setrecursionlimit(10000)


def create_blank_morphology_dict():
    atom_dictionary = {
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
        "xy": 0.0,
        "xz": 0.0,
        "yz": 0.0,
        "lx": 0.0,
        "ly": 0.0,
        "lz": 0.0,
        "natoms": 0,
    }
    return atom_dictionary


def find_index(string, character):
    """
    This function returns the locations of an inputted character in an inputted string
    """
    index = 0
    locations = []
    while index < len(string):
        if string[index] == character:
            locations.append(index)
        index += 1
    if len(locations) == 0:
        return None
    return locations


def calculate_separation(atom1, atom2):
    """
    This function calculates the distance between two input points (either as lists or
    np.arrays)
    """
    atom1 = np.array(atom1)
    atom2 = np.array(atom2)
    return np.sqrt(np.sum((atom1 - atom2) ** 2))


def calc_COM(list_of_positions, list_of_atom_types=None, list_of_atom_masses=None):
    """
    This function calculates the centre of mass of a collection of sites/atoms
    (listOfPositions) with corresponding type (listOfAtomTypes) or mass (listOfMasses)
    if list_of_atom_masses is not specified, then list_of_atom_types must be.
    """
    mass_weighted = np.array([0.0, 0.0, 0.0])
    if list_of_atom_masses is None:
        list_of_atom_masses = []
        for atom_type in list_of_atom_types:
            # Masses obtained from nist.gov, for the atoms we are likely to
            # simulate the most.
            # Add in new atoms here if your molecule requires it!
            if atom_type.lower()[:2] == "si":
                list_of_atom_masses.append(27.976926)
            elif atom_type.lower()[:2] == "mo":
                list_of_atom_masses.append(95.960000)
            elif atom_type.lower()[:2] == "nb":
                list_of_atom_masses.append(92.906380)
            elif atom_type.lower()[:2] == "te":
                list_of_atom_masses.append(127.60000)
            elif atom_type.lower()[:2] == "ni":
                list_of_atom_masses.append(140.91120)
            elif atom_type.lower()[:2] == "ga":
                list_of_atom_masses.append(69.723000)
            elif atom_type.lower()[:2] == "mn":
                list_of_atom_masses.append(54.938045)
            elif atom_type.lower()[:2] == "cu":
                list_of_atom_masses.append(63.546000)
            elif atom_type.lower()[:2] == "ag":
                list_of_atom_masses.append(107.86820)
            elif atom_type.lower()[:2] == "au":
                list_of_atom_masses.append(196.96657)
            elif atom_type.lower()[0] == "c":
                list_of_atom_masses.append(12.000000)
            elif atom_type.lower()[0] == "h":
                list_of_atom_masses.append(1.0078250)
            elif atom_type.lower()[0] == "s":
                list_of_atom_masses.append(31.972071)
            elif atom_type.lower()[0] == "o":
                list_of_atom_masses.append(15.994914)
            elif atom_type.lower()[0] == "n":
                list_of_atom_masses.append(14.003074)
            elif atom_type.lower()[0] == "v":
                list_of_atom_masses.append(50.941500)
            elif atom_type.lower()[0] == "f":
                list_of_atom_masses.append(18.998403)
            else:
                print("Unknown atomic mass {:s}. Setting as 1.0.".format(atom_type))
                list_of_atom_masses.append(1.0)
    total_mass = np.sum(list_of_atom_masses)
    for atom_ID, position in enumerate(list_of_positions):
        for axis in range(3):
            mass_weighted[axis] += position[axis] * list_of_atom_masses[atom_ID]
    return mass_weighted / float(total_mass)


def find_axis(atom1, atom2, normalise=True):
    """
    This function determines the normalised vector from the location of atom1 to atom2.
    The positions can enter as lists or arrays, but are output as arrays
    """
    x_sep = atom2[0] - atom1[0]
    y_sep = atom2[1] - atom1[1]
    z_sep = atom2[2] - atom1[2]
    if normalise is True:
        axis_vector = normalise_vec(np.array([x_sep, y_sep, z_sep]))
    else:
        axis_vector = np.array([x_sep, y_sep, z_sep])
    return axis_vector


def normalise_vec(vector):
    """
    This function normalises an input vector to unit magnitude
    """
    return vector / float(np.sqrt(np.sum(vector) ** 2))


def get_rotation_matrix(vector1, vector2):
    """
    This function returns the rotation matrix around the origin that maps vector1 to
    vector 2
    """
    cross_product = np.cross(vector1, vector2)
    sin_angle = np.sqrt(
        (
            (cross_product[0] ** 2)
            + ((cross_product[1]) ** 2)
            + ((cross_product[2]) ** 2)
        )
    )
    cos_angle = np.dot(vector1, vector2)
    skew_matrix = np.array(
        [
            [0, -cross_product[2], cross_product[1]],
            [cross_product[2], 0, -cross_product[0]],
            [-cross_product[1], cross_product[0], 0],
        ]
    )
    skew_matrix_squared = skew_matrix @ skew_matrix
    rot_matrix = (
        np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        + skew_matrix
        + skew_matrix_squared * ((1 - cos_angle) / (sin_angle ** 2))
    )
    return rot_matrix


def parallel_sort(list1, list2):
    """
    This function sorts a pair of lists by the first list in ascending order (for
    example, atom mass and corresponding position can be input, sorted by ascending
    mass, and the two lists output, where the mass[atom_i] still corresponds to
    position[atom_i]
    """
    list1, list2 = zip(*sorted(zip(list1, list2)))
    return list1, list2


def write_CSV(file_name, data):
    """
    Writes a CSV file given a 2D array `data' of arbitrary size
    """
    with open(file_name, "w+") as csv_file:
        document = csv.writer(csv_file, delimiter=",")
        for row in data:
            document.writerow(list(row))
    print("CSV written to", file_name)


def add_unwrapped_positions(input_dictionary):
    """
    This function takes a run_HOOMD.py input dictionary and updates the
    'unwrapped_position' key based on the values of the 'position' and 'image' keys
    """
    simulation_dimensions = [
        input_dictionary["lx"],
        input_dictionary["ly"],
        input_dictionary["lz"],
    ]
    input_dictionary["unwrapped_position"] = [0] * len(input_dictionary["position"])
    for i in range(len(input_dictionary["position"])):
        position = input_dictionary["position"][i]
        if len(input_dictionary["image"]) > 0:
            image = input_dictionary["image"][i]
        else:
            image = [0, 0, 0]
        unwrapped_position = []
        for axis in range(len(image)):
            unwrapped_position.append(
                (image[axis] * simulation_dimensions[axis]) + position[axis]
            )
        input_dictionary["unwrapped_position"][i] = unwrapped_position
    return input_dictionary


def replace_wrapped_positions(input_dictionary):
    """
    This function takes a MorphCT input dictionary and replaces the 'position' and
    'image' keys with the 'unwrapped_position' key and '[0, 0, 0]' respectively.
    """
    for atom_ID, unwrapped_position in enumerate(
        input_dictionary["unwrapped_position"]
    ):
        input_dictionary["position"][atom_ID] = unwrapped_position
        input_dictionary["image"][atom_ID] = [0, 0, 0]
    return input_dictionary


def add_wrapped_positions(input_dictionary):
    """
    This function takes a run_HOOMD.py input dictionary and updates the 'position' and
    'image' keys based on the values of the 'unwrapped_position' key
    """
    simulation_dimensions = [
        input_dictionary["lx"],
        input_dictionary["ly"],
        input_dictionary["lz"],
    ]
    input_dictionary["position"] = [0] * len(input_dictionary["unwrapped_position"])
    input_dictionary["image"] = [0] * len(input_dictionary["unwrapped_position"])
    for atom_ID in range(len(input_dictionary["unwrapped_position"])):
        position = copy.deepcopy(input_dictionary["unwrapped_position"][atom_ID])
        image_coords = [0, 0, 0]
        for axis in range(len(position)):
            if position[axis] > (simulation_dimensions[axis] / 2.0):
                while position[axis] > (simulation_dimensions[axis] / 2.0):
                    image_coords[axis] += 1
                    position[axis] -= simulation_dimensions[axis]
            elif position[axis] < -(simulation_dimensions[axis] / 2.0):
                while position[axis] < -(simulation_dimensions[axis] / 2.0):
                    image_coords[axis] -= 1
                    position[axis] += simulation_dimensions[axis]
        input_dictionary["position"][atom_ID] = position
        input_dictionary["image"][atom_ID] = image_coords
    return input_dictionary


def add_masses(input_dictionary):
    """
    This function takes a run_HOOMD.py input dictionary and updates the 'mass' key based
    on the values of the 'type' key. Note that more hardcoding is required to add
    aditional atom types
    """
    input_dictionary["mass"] = [1.0] * len(input_dictionary["type"])
    for atom_ID in range(len(input_dictionary["type"])):
        if "H" in input_dictionary["type"][atom_ID]:
            input_dictionary["mass"][atom_ID] = 1.007825
        elif "C" in input_dictionary["type"][atom_ID]:
            input_dictionary["mass"][atom_ID] = 12.00000
        elif "S" in input_dictionary["type"][atom_ID]:
            input_dictionary["mass"][atom_ID] = 31.972071
        elif "O" in input_dictionary["type"][atom_ID]:
            input_dictionary["mass"][atom_ID] = 15.994914
        elif "N" in input_dictionary["type"][atom_ID]:
            input_dictionary["mass"][atom_ID] = 14.003074
    return input_dictionary


def add_diameters(input_dictionary):
    """
    This function takes a run_HOOMD.py input dictionary and updates the 'diameter' key
    based on the values of the 'type' key. Values are given in A. Note that more
    hardcoding is required to add aditional atom types
    """
    input_dictionary["diameter"] = [1.0] * len(input_dictionary["type"])
    for atom_ID in range(len(input_dictionary["type"])):
        if "H" in input_dictionary["type"][atom_ID]:
            input_dictionary["diameter"][atom_ID] = 1.200
        elif "C" in input_dictionary["type"][atom_ID]:
            input_dictionary["diameter"][atom_ID] = 1.700
        elif "S" in input_dictionary["type"][atom_ID]:
            input_dictionary["diameter"][atom_ID] = 1.800
        elif "O" in input_dictionary["type"][atom_ID]:
            input_dictionary["diameter"][atom_ID] = 1.520
        elif "N" in input_dictionary["type"][atom_ID]:
            input_dictionary["diameter"][atom_ID] = 1.550
    return input_dictionary


def get_terminating_positions(
    current_atom_posn, bonded_atom_positions, number_of_units_to_add
):
    # Given a current_atom_posn and several bonded_atom_positions we can add
    # number_of_units_to_add different terminating units to the current_atom
    # through a series of geometric checks. First get the vector to the average
    # position of the bonded neighbours
    hydrogen_positions = []
    average_position_of_bonded_atoms = np.array([0.0, 0.0, 0.0])
    for bonded_atom_posn in bonded_atom_positions:
        bond_vector = np.array(bonded_atom_posn) - current_atom_posn
        bond_vector /= np.linalg.norm(bond_vector)
        average_position_of_bonded_atoms += bond_vector
    [x, y, z] = current_atom_posn + (
        -1.06
        * (
            average_position_of_bonded_atoms
            / np.linalg.norm(average_position_of_bonded_atoms)
        )
    )
    if number_of_units_to_add == 1:
        # Easy, this is the perylene code
        # Simply reverse the bonded vector and make it the hydrogen position at
        # a distance of 1.06 angstroems
        hydrogen_positions.append(np.array([x, y, z]))
    # Initial position for all hydrogens
    elif number_of_units_to_add == 2:
        # As above (to get the right plane), but then rotated +(109.5/2) degrees
        # and -(109.5/2) degrees around the bonding axis
        rotation_axis = np.array(bonded_atom_positions[0]) - np.array(
            bonded_atom_positions[-1]
        )
        rotation_axis /= np.linalg.norm(rotation_axis)
        # Rotation matrix calculations from:
        # http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
        # The array that describes the 3D rotation of (x, y, z) around the point
        # (a, b, c) through
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
            hydrogen_positions.append(new_position)
    elif number_of_units_to_add == 3:
        # As for one (to get the right side of the bonded atom), rotate the
        # first one up by 70.5 (180 - 109.5) and then rotate around by 109.5
        # degrees for the other two. The first hydrogen can be rotated around
        # any axis perpendicular to the only bond present
        axis_to_bond = current_atom_posn - np.array(bonded_atom_positions[0])
        # Now find one of the set of vectors [i, j, k] perpendicular to this one
        # so we can place the first hydrogen.
        # Do this by setting i = j = 1 and solve for k (given that
        # current_atom_posn[0]*i + current_atom_posn[1]*j
        # + current_atom_posn[2]*k = 0)
        first_hydrogen_rotation_axis = np.array(
            [1, 1, -(axis_to_bond[0] + axis_to_bond[1]) / axis_to_bond[2]]
        )
        first_hydrogen_rotation_axis /= np.linalg.norm(first_hydrogen_rotation_axis)

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
        hydrogen_positions.append(new_position)
        # Second and third hydrogens
        # Rotate these from the new_position +/-120 degrees around the vector
        # axis_to_bond from the position current_atom_posn - axis_to_bond
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
            hydrogen_positions.append(new_hydrogen_position)
    return hydrogen_positions


def load_morphology_xml(xml_path, sigma=1.0):
    # xml has SimDims as <box
    # Positions as <position and <image
    # Velocities as <velocity
    # Mass as <mass
    # Diameters as <diameter
    # Type as <type
    # "Body" as <body
    # Bonds as <bond, with each line as bondA, bondB, etc.
    # Angles as <angle, with each angle as angleA, angleB, etc.
    # Dihedral as <dihedral
    # Improper as <improper (usually none in xml)
    # Charge as <charge
    # Set default tilts so when we use simdims later they exisit always
    atom_dictionary = create_blank_morphology_dict()
    record = False
    with open(xml_path, "r") as xml_file:
        xml_data = xml_file.readlines()
        for line in xml_data:
            if ("</" in line) or ("<!--" in line):
                record = False
            elif ("<configuration" in line) or ("<box" in line):
                # Get configuration data from this line (timestep, natoms etc)
                split_line = line.split(" ")
                for i in range(1, len(split_line)):
                    equals_loc = find_index(split_line[i], "=")
                    if (equals_loc is None) or ("units" in split_line[i]):
                        # Skip any elements without equals or to do with the
                        # distance units
                        continue
                    quot_loc = find_index(split_line[i], '"')
                    if "." in split_line[i][quot_loc[0] + 1 : quot_loc[1]]:
                        # Catch float in the value (excludes the = and quotation
                        # marks)
                        atom_dictionary[split_line[i][: equals_loc[0]].lower()] = float(
                            split_line[i][quot_loc[0] + 1 : quot_loc[1]]
                        )
                    else:
                        atom_dictionary[split_line[i][: equals_loc[0]].lower()] = int(
                            split_line[i][quot_loc[0] + 1 : quot_loc[1]]
                        )
            elif "<position" in line:
                record = True
                record_type = "position"
                continue
            elif "<image" in line:
                record = True
                record_type = "image"
                continue
            elif "<mass" in line:
                record = True
                record_type = "mass"
                continue
            elif "<diameter" in line:
                record = True
                record_type = "diameter"
                continue
            elif "<type" in line:
                record = True
                record_type = "type"
                continue
            elif "<body" in line:
                record = True
                record_type = "body"
                continue
            elif ("<bond" in line) and ("_coeff" not in line):
                record = True
                record_type = "bond"
                continue
            elif ("<angle" in line) and ("_coeff" not in line):
                record = True
                record_type = "angle"
                continue
            elif ("<dihedral" in line) and ("_coeff" not in line):
                record = True
                record_type = "dihedral"
                continue
            elif ("<improper" in line) and ("_coeff" not in line):
                record = True
                record_type = "improper"
                continue
            elif "<charge" in line:
                record = True
                record_type = "charge"
                continue
            # Now we know what the variable is, append it to the dictionary data
            if record is True:
                split_line = line.split()
                if record_type == "position":
                    # NOTE: VELOCITIES ARE NOT NORMALISED IN THE MORPHOLOGY
                    # FILE...DO THEY NEED TO BE SCALED BY SIGMA OR NOT?
                    # CURRENTLY THEY ARE.
                    # Write to dictionary as floats scaled by sigma
                    if len(split_line) == 1:
                        atom_dictionary[record_type].append(float(split_line[0]))
                        continue
                    for i in range(len(split_line)):
                        split_line[i] = float(split_line[i])
                    atom_dictionary[record_type].append(split_line)
                elif (
                    (record_type == "mass")
                    or (record_type == "diameter")
                    or (record_type == "charge")
                ):
                    # Write to dictionary as floats
                    if len(split_line) == 1:
                        atom_dictionary[record_type].append(float(split_line[0]))
                        continue
                    for i in range(len(split_line)):
                        split_line[i] = float(split_line[i])
                    atom_dictionary[record_type].append(split_line)
                elif (record_type == "image") or (record_type == "body"):
                    # Write to dictionary as int
                    if len(split_line) == 1:
                        atom_dictionary[record_type].append(int(split_line[0]))
                        continue
                    for i in range(len(split_line)):
                        split_line[i] = int(split_line[i])
                    atom_dictionary[record_type].append(split_line)
                elif record_type == "type":
                    # Write to dictionary as str
                    atom_dictionary[record_type].append(str(split_line[0]))
                else:
                    # (recordType == 'bond') or (recordType == 'angle') or
                    # (recordType == 'dihedral') or (recordType == 'improper')
                    # Write to dictionary as combination
                    split_line[0] = str(split_line[0])
                    for i in range(1, len(split_line)):
                        split_line[i] = int(split_line[i])
                    atom_dictionary[record_type].append(split_line)
    if sigma != 1.0:
        atom_dictionary = scale(atom_dictionary, sigma)
    return atom_dictionary


def load_FF_xml(xml_path, mapping=False):
    FF_dict = {
        "lj": [],
        "dpd": [],
        "bond": [],
        "angle": [],
        "dihedral": [],
        "improper": [],
    }
    with open(xml_path, "r") as xml_file:
        xml_data = xml_file.readlines()
        for line in xml_data:
            if "</" in line:
                record = False
            elif "<lj" in line:
                record = True
                record_type = "lj"
                continue
            elif "<dpd" in line:
                record = True
                record_type = "dpd"
                continue
            elif "<bond" in line:
                record = True
                record_type = "bond"
                continue
            elif "<angle" in line:
                record = True
                record_type = "angle"
                continue
            elif "<dihedral" in line:
                record = True
                record_type = "dihedral"
                continue
            elif "<improper" in line:
                record = True
                record_type = "improper"
                continue
            # Now we know what the variable is, append it to the dictionary data
            if record is True:
                # Write to dictionary as combination
                split_line = line.split(" ")
                # Remove the "\n"
                split_line[-1] = split_line[-1][:-1]
                split_line[0] = str(split_line[0])
                for i in range(1, len(split_line)):
                    split_line[i] = float(split_line[i])
                FF_dict[record_type].append(split_line)
    # Now remap the names of the constraints if any mappings have been specified
    if mapping is not False:
        for constraint_type in list(FF_dict.keys()):
            for index, constraint in enumerate(FF_dict[constraint_type]):
                # Split the constraint name up based on each atom type
                constraint_name = copy.deepcopy(constraint[0].split("-"))
                # Remap each atom
                for atom_loc, atom_type in enumerate(constraint_name):
                    constraint_name[atom_loc] = mapping[atom_type]
                # Apply the mapping to the FFDict
                FF_dict[constraint_type][index][0] = "-".join(constraint_name)
    return FF_dict


def check_constraint_names(AA_morphology_dict):
    # A function that renames the constraints based on the atom types given in
    # the dictionary
    constraint_types = ["bond", "angle", "dihedral", "improper"]
    for constraint_type in constraint_types:
        for constraint_ID, constraint in enumerate(AA_morphology_dict[constraint_type]):
            new_constraint_name = "-".join(
                [AA_morphology_dict["type"][atom_ID] for atom_ID in constraint[1:]]
            )
            # Update the dict if the name has changed
            if constraint[0] != new_constraint_name:
                AA_morphology_dict[constraint_type][constraint_ID][
                    0
                ] = new_constraint_name
    return AA_morphology_dict


def write_morphology_xml(
    input_dictionary, output_file, sigma=1.0, check_wrapped_posns=True
):
    # Firstly, scale everything by the inverse of the provided sigma value
    box_lengths = ["lx", "ly", "lz"]
    tilt_factors = ["xy", "yz", "xz"]
    if sigma != 1.0:
        input_dictionary = scale(input_dictionary, 1.0 / sigma)
    # Now need to check the positions of the atoms to ensure that everything is
    # correctly contained inside the box
    if check_wrapped_posns is True:
        if all(
            np.isclose(
                [input_dictionary[box_length] for box_length in box_lengths], [0, 0, 0]
            )
        ):
            print(
                "No box length specified, cannot wrap positions. Continuing as if"
                " check_wrapped_posns is False..."
            )
        else:
            if any(
                [tilt_factor in input_dictionary.keys() for tilt_factor in tilt_factors]
            ) and any(
                [input_dictionary[tilt_factor] != 0 for tilt_factor in tilt_factors]
            ):
                print("Can't check atom wrapping for cells with a non-zero tilt factor")
            else:
                print("Checking wrapped positions before writing xml...")
                input_dictionary = check_wrapped_positions(input_dictionary)
    # Need to add natoms if it doesn't exist already
    if "natoms" not in input_dictionary.keys():
        input_dictionary["natoms"] = len(input_dictionary["position"])
    # Add Boiler Plate first
    lines_to_write = [
        '<?xml version="1.0" encoding="UTF-8"?>\n',
        '<hoomd_xml version="1.4">\n',
        '<configuration time_step="0" dimensions="3" natoms="{:d}" >\n'.format(
            input_dictionary["natoms"]
        ),
        '<box lx="{0:f}" ly="{1:f}" lz="{2:f}"'.format(
            input_dictionary["lx"], input_dictionary["ly"], input_dictionary["lz"]
        ),
    ]
    if all([tilt_factor in input_dictionary.keys() for tilt_factor in tilt_factors]):
        lines_to_write[-1] = "".join(
            [
                lines_to_write[-1],
                ' xy="{0:f}" xz="{1:f}" yz="{2:f}" />\n'.format(
                    input_dictionary["xy"],
                    input_dictionary["xz"],
                    input_dictionary["yz"],
                ),
            ]
        )
    else:
        lines_to_write[-1] = "".join([lines_to_write[-1], '" />\n'])
    # Position
    lines_to_write.append('<position num="{:d}">\n'.format(input_dictionary["natoms"]))
    for position_data in input_dictionary["position"]:
        lines_to_write.append(" ".join(str(coord) for coord in position_data) + "\n")
    lines_to_write.append("</position>\n")
    # Image
    lines_to_write.append('<image num="{:d}">\n'.format(input_dictionary["natoms"]))
    for image_data in input_dictionary["image"]:
        lines_to_write.append(" ".join(str(coord) for coord in image_data) + "\n")
    lines_to_write.append("</image>\n")
    # Mass
    lines_to_write.append('<mass num="{:d}">\n'.format(input_dictionary["natoms"]))
    for mass_data in input_dictionary["mass"]:
        lines_to_write.append("".join([str(mass_data), "\n"]))
    lines_to_write.append("</mass>\n")
    # Diameter
    lines_to_write.append('<diameter num="{:d}">\n'.format(input_dictionary["natoms"]))
    for diameter_data in input_dictionary["diameter"]:
        lines_to_write.append("".join([str(diameter_data), "\n"]))
    lines_to_write.append("</diameter>\n")
    # Type
    lines_to_write.append('<type num="{:d}">\n'.format(input_dictionary["natoms"]))
    for type_data in input_dictionary["type"]:
        lines_to_write.append("".join([str(type_data), "\n"]))
    lines_to_write.append("</type>\n")
    # Body
    lines_to_write.append('<body num="{:d}">\n'.format(input_dictionary["natoms"]))
    for body_data in input_dictionary["body"]:
        lines_to_write.append("".join([str(body_data), "\n"]))
    lines_to_write.append("</body>\n")
    # Bond
    lines_to_write.append('<bond num="{:d}">\n'.format(len(input_dictionary["bond"])))
    for bond_data in input_dictionary["bond"]:
        lines_to_write.append(" ".join(str(coord) for coord in bond_data) + "\n")
    lines_to_write.append("</bond>\n")
    # Angle
    lines_to_write.append('<angle num="{:d}">\n'.format(len(input_dictionary["angle"])))
    for angle_data in input_dictionary["angle"]:
        lines_to_write.append(" ".join(str(coord) for coord in angle_data) + "\n")
    lines_to_write.append("</angle>\n")
    # Dihedral
    lines_to_write.append(
        '<dihedral num="{:d}">\n'.format(len(input_dictionary["dihedral"]))
    )
    for dihedral_data in input_dictionary["dihedral"]:
        lines_to_write.append(" ".join(str(coord) for coord in dihedral_data) + "\n")
    lines_to_write.append("</dihedral>\n")
    # Improper
    lines_to_write.append(
        '<improper num="{:d}">\n'.format(len(input_dictionary["improper"]))
    )
    for improper_data in input_dictionary["improper"]:
        lines_to_write.append(" ".join(str(coord) for coord in improper_data) + "\n")
    lines_to_write.append("</improper>\n")
    # Charge
    lines_to_write.append('<charge num="{:d}">\n'.format(input_dictionary["natoms"]))
    for charge_data in input_dictionary["charge"]:
        lines_to_write.append("".join([str(charge_data), "\n"]))
    lines_to_write.append("</charge>\n")
    lines_to_write.append("</configuration>\n")
    lines_to_write.append("</hoomd_xml>\n")
    with open(output_file, "w+") as xml_file:
        xml_file.writelines(lines_to_write)
    print("XML file written to", str(output_file))


def write_xyz_file(input_dict, output_file):
    """
    This function takes an input dictionary and converts it to an xyz for use in DFT
    calculations
    """
    # First line is atom numbers, second line is boiler plate
    rows_to_write = [
        "{:d}\n".format(input_dict["natoms"]),
        "xyz file generated from xml using helper_functions.xml_to_xyz\n",
    ]
    # Format of xyz is Type, X Pos, Y Pos, Z Pos
    for atom_ID in range(len(input_dict["type"])):
        # Note, this will break for atoms that have two-letter symbols (e.g. Al,
        # Ca etc.)
        atom_type = input_dict["type"][atom_ID][0]
        while len(atom_type) < 10:
            atom_type += " "
        atom_X = str(input_dict["position"][atom_ID][0])
        while len(atom_X) < 20:
            atom_X += " "
        atom_Y = str(input_dict["position"][atom_ID][1])
        while len(atom_Y) < 20:
            atom_Y += " "
        atom_Z = str(input_dict["position"][atom_ID][2])
        line_to_write = "".join([atom_type, atom_X, atom_Y, atom_Z, "\n"])
        rows_to_write.append(line_to_write)
    with open(output_file, "w+") as xyz_file:
        xyz_file.writelines(rows_to_write)
    print("XYZ data written to", str(output_file))


def increment_atom_IDs(
    original_input_dictionary,
    ghost_dictionary,
    increment,
    modify_ghost_dictionary=False,
):
    input_dictionary = copy.deepcopy(original_input_dictionary)
    con_types = ["bond", "angle", "dihedral", "improper"]
    for con_type in con_types:
        for constraint_no, constraint in enumerate(input_dictionary[con_type]):
            input_dictionary[con_type][constraint_no][1:] = [
                x + increment for x in input_dictionary[con_type][constraint_no][1:]
            ]
    if modify_ghost_dictionary is True:
        for bond_no, bond in enumerate(ghost_dictionary["bond"]):
            if str(bond[1])[0] == "_":
                ghost_dictionary["bond"][bond_no][1] = int(bond[1][1:]) + increment
            if str(bond[2])[0] == "_":
                ghost_dictionary["bond"][bond_no][2] = int(bond[2][1:]) + increment
    return input_dictionary, ghost_dictionary


def scale(input_dictionary, scale_factor):
    for ID, position in enumerate(input_dictionary["position"]):
        input_dictionary["position"][ID] = list(scale_factor * np.array(position))
    for element in ["lx", "ly", "lz"]:
        if element in input_dictionary:
            input_dictionary[element] *= scale_factor
    return input_dictionary


def rotate(
    input_dictionary, theta, rotate_around_point=[0, 0, 0], rotate_around_axis=[0, 0, 1]
):
    input_dictionary = add_unwrapped_positions(input_dictionary)
    rotate_around_axis = list(
        np.array(rotate_around_axis) / np.linalg.norm(rotate_around_axis)
    )
    # Rotation matrix calculations from:
    # http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
    # The array that describes the 3D rotation of (x, y, z) around the point
    # (a, b, c) through the unit axis <u, v, w> by the angle theta is given by:
    # [ (a(v^2 + w^2) - u(bv + cw - ux - vy - wz))(1 - cos(theta))
    #    + x*cos(theta) + (-cv + bw - wy + vz)sin(theta),
    #   (b(u^2 + w^2) - v(au + cw - ux - vy - wz))(1 - cos(theta))
    #    + y*cos(theta) + (cu - aw + wx - uz)sin(theta),
    #   (c(u^2 + v^2) - w(au + bv - ux - vy - wz))(1 - cos(theta))
    #    + z*cos(theta) + (-bu + av - vx + uy)sin(theta) ]
    # DEFAULT BEHAVIOUR: Rotate the entire dictionary by theta around the z-axis
    # centred at the origin
    for AAID, [x, y, z] in enumerate(input_dictionary["unwrapped_position"]):
        [a, b, c] = rotate_around_point
        [u, v, w] = rotate_around_axis
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
        input_dictionary["unwrapped_position"][AAID] = list(new_position)
    # All the images are probably messed up now, so fix that
    input_dictionary = replace_wrapped_positions(input_dictionary)
    return input_dictionary


def centre(input_dictionary, centre_of_mass):
    COM = np.array(centre_of_mass)
    for index, position in enumerate(input_dictionary["position"]):
        input_dictionary["position"][index] = list(position - COM)
    return input_dictionary


def check_wrapped_positions(input_dictionary):
    atom_positions = np.array(input_dictionary["position"])
    atom_images = np.array(input_dictionary["image"])
    xhi = input_dictionary["lx"] / 2.0
    xlo = -input_dictionary["lx"] / 2.0
    yhi = input_dictionary["ly"] / 2.0
    ylo = -input_dictionary["ly"] / 2.0
    zhi = input_dictionary["lz"] / 2.0
    zlo = -input_dictionary["lz"] / 2.0
    for atom_ID in range(len(atom_positions)):
        while atom_positions[atom_ID][0] > xhi:
            atom_positions[atom_ID][0] -= input_dictionary["lx"]
            atom_images[atom_ID][0] += 1
        while atom_positions[atom_ID][0] < xlo:
            atom_positions[atom_ID][0] += input_dictionary["lx"]
            atom_images[atom_ID][0] -= 1
        while atom_positions[atom_ID][1] > yhi:
            atom_positions[atom_ID][1] -= input_dictionary["ly"]
            atom_images[atom_ID][1] += 1
        while atom_positions[atom_ID][1] < ylo:
            atom_positions[atom_ID][1] += input_dictionary["ly"]
            atom_images[atom_ID][1] -= 1
        while atom_positions[atom_ID][2] > zhi:
            atom_positions[atom_ID][2] -= input_dictionary["lz"]
            atom_images[atom_ID][2] += 1
        while atom_positions[atom_ID][2] < zlo:
            atom_positions[atom_ID][2] += input_dictionary["lz"]
            atom_images[atom_ID][2] -= 1
    input_dictionary["position"] = list(atom_positions)
    input_dictionary["image"] = list(atom_images)
    return input_dictionary


def get_CPU_cores():
    # Determine the number of available processors, either by querying the
    # SLURM_NPROCS environment variable, or by using multiprocessing to count
    # the number of visible CPUs.
    try:
        proc_IDs = list(np.arange(int(os.environ.get("SLURM_NPROCS"))))
    except (AttributeError, TypeError):
        # Was not loaded using SLURM, so use all physical processors
        proc_IDs = list(np.arange(mp.cpu_count()))
    return proc_IDs


def write_to_file(log_file, string_list, mode="log_file"):
    if mode == "output_file":
        open_as = "w+"
    else:
        open_as = "a+"
    if log_file == "stdout":
        if sys.stdout is not None:
            for line in string_list:
                sys.stdout.writelines("".join([line, "\n"]))
    else:
        with open(log_file, open_as) as log_write:
            for line in string_list:
                log_write.writelines("".join([line, "\n"]))


def load_pickle(pickle_location):
    print("".join(["Loading Pickle from ", str(pickle_location), "..."]))
    try:
        with open(pickle_location, "rb") as pickle_file:
            objects = pickle.load(pickle_file)
    except UnicodeDecodeError:  # Python 2/3 fix
        print("Old pickle! Loading it using Python 2...")
        with open(pickle_location, "rb") as pickle_file:
            objects = pickle.load(pickle_file, encoding="latin1")
    print("Pickle loaded successfully!")
    return [objects[0], objects[1], objects[2], objects[3], objects[4]]


def write_pickle(to_pickle, pickle_file_name):
    print("Writing pickle file...")
    with open(pickle_file_name, "wb+") as pickle_file:
        pickle.dump(to_pickle, pickle_file)
    print("Pickle file written to", pickle_file_name)


def obtain_bonded_list(bond_list):
    # Create a lookup table `neighbour list' for all connected atoms called
    # {bondedAtoms}
    bonded_atoms = {}
    for bond in bond_list:
        if bond[1] not in bonded_atoms:
            bonded_atoms[bond[1]] = [bond[2]]
        else:
            bonded_atoms[bond[1]].append(bond[2])
        if bond[2] not in bonded_atoms:
            bonded_atoms[bond[2]] = [bond[1]]
        else:
            bonded_atoms[bond[2]].append(bond[1])
    return bonded_atoms


def convert_string_to_int(x):
    for i in range(len(x)):
        try:
            return int(x[i:])
        except:
            continue
    return 99999


def fix_images(original_morphology):
    def check_bonds(morphology, bond_dict):
        for bond in morphology["bond"]:
            posn1 = np.array(morphology["position"][bond[1]]) + (
                np.array(morphology["image"][bond[1]])
                * np.array([morphology["lx"], morphology["ly"], morphology["lz"]])
            )
            posn2 = np.array(morphology["position"][bond[2]]) + (
                np.array(morphology["image"][bond[2]])
                * np.array([morphology["lx"], morphology["ly"], morphology["lz"]])
            )
            separation = calculate_separation(posn1, posn2)
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
        for atom_ID, image in enumerate(morphology["image"]):
            if image != [0, 0, 0]:
                morphology["image"][atom_ID] = [0, 0, 0]
        return morphology

    def get_bond_dict(morphology):
        bond_dict = {
            atom_ID: [] for atom_ID, atom_type in enumerate(morphology["type"])
        }
        for bond in morphology["bond"]:
            bond_dict[bond[1]].append(bond[2])
            bond_dict[bond[2]].append(bond[1])
        return bond_dict

    def move_bonded_atoms(central_atom, morphology, bond_dict):
        for bonded_atom in bond_dict[central_atom]:
            atom1posn = morphology["position"][central_atom]
            atom2posn = morphology["position"][bonded_atom]
            sep_vec = np.array(atom1posn) - np.array(atom2posn)
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
            else:
                pass
        return morphology

    zeroed_morphology = zero_out_images(original_morphology)
    bond_dict = get_bond_dict(zeroed_morphology)
    fixed_morphology = check_bonds(zeroed_morphology, bond_dict)
    return fixed_morphology


# ---============================---
# ---=== KMC HELPER FUNCTIONS ===---
# ---============================---
def calculate_carrier_hop_rate(
    lambda_ij,
    T_ij,
    delta_E_ij,
    prefactor,
    temp,
    use_VRH=False,
    rij=0.0,
    VRH_delocalisation=1.0,
    boltz_pen=False,
):
    # Based on the input parameters, can make this the semiclassical Marcus
    # Hopping Rate Equation, or a more generic Miller Abrahams-based hop
    # Firstly, to prevent divide-by-zero errors:
    if T_ij == 0.0:
        return 0
    # Regardless of hopping type, sort out the prefactor first:
    k_ij = (
        prefactor
        * ((2 * np.pi) / hbar)
        * (T_ij ** 2)
        * np.sqrt(1.0 / (4 * lambda_ij * np.pi * k_B * temp))
    )
    # VRH?
    if use_VRH is True:
        k_ij *= np.exp(-(rij / VRH_delocalisation))
    # Simple Boltzmann energy penalty?
    if boltz_pen is True:
        # Only apply the penalty if delta_E_ij is positive
        if delta_E_ij > 0.0:
            k_ij *= np.exp(-(delta_E_ij / (k_B * temp)))
        # Otherwise, k_ij *= 1
    else:
        k_ij *= np.exp(-((delta_E_ij + lambda_ij) ** 2) / (4 * lambda_ij * k_B * temp))
    return k_ij


def calculate_FRET_hop_rate(prefactor, lifetime_parameter, r_F, rij, delta_E_ij, T):
    # Foerster Transport Hopping Rate Equation
    # The prefactor included here is a bit of a bodge to try and get the
    # mean-free paths of the excitons more in line with the 5nm of experiment.
    # Possible citation: 10.3390/ijms131217019 (they do not do the simulation
    # they just point out some limitations of FRET which assumes point-dipoles
    # which does not necessarily work in all cases)
    if delta_E_ij <= 0:
        boltzmann_factor = 1
    else:
        boltzmann_factor = np.exp(-(elementary_charge * delta_E_ij) / (k_B * T))
    k_FRET = prefactor * (1 / lifetime_parameter) * (r_F / rij) ** 6 * boltzmann_factor
    return k_FRET


def calculate_miller_abrahams_hop_rate(prefactor, separation, radius, delta_E_ij, T):
    k_ij = prefactor * np.exp(-2 * separation / radius)
    if delta_E_ij > 0:
        k_ij *= np.exp(-delta_E_ij / (k_B * T))
    return k_ij


def determine_event_tau(
    rate,
    event_type="None",
    slowest_event=None,
    fastest_event=None,
    maximum_attempts=None,
    log_file=None,
):
    # Use the KMC algorithm to determine the wait time to this hop
    if rate != 0:
        counter = 0
        while True:
            if maximum_attempts is not None:
                # Write an error if we've hit the maximum number of attempts
                if counter == maximum_attempts:
                    if "hop" in event_type:
                        return None
                    else:
                        if ("injection" not in event_type) and (log_file is not None):
                            write_to_file(
                                log_file,
                                [
                                    "Attempted {0:d} times to obtain a {1:s}-type event timescale withing the tolerances: {2:.2e} <= tau < {3:.2e} with the given rate {4:.2e}, all without success. Permitting the event anyway with the next random number...".format(
                                        maximum_attempts,
                                        event_type,
                                        fastest_event,
                                        slowest_event,
                                        rate,
                                    )
                                ],
                            )

            x = np.random.random()
            # Ensure that we don't get exactly 0.0 or 1.0, which would break our
            # logarithm
            if (x == 0.0) or (x == 1.0):
                continue
            tau = -np.log(x) / rate
            if (
                (fastest_event is not None)
                and (slowest_event is not None)
                and (maximum_attempts is not None)
            ):
                if ((tau > fastest_event) and (tau < slowest_event)) or (
                    counter == maximum_attempts
                ):
                    break
                else:
                    counter += 1
                    continue
            break
    else:
        # If rate == 0, then make the hopping time extremely long
        tau = 1E99
    return tau
