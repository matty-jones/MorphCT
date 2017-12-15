import sys
import numpy as np
sys.path.append('../../code')
from morphct.code import helperFunctions


def check_bonds(morphology, bond_dict):
    periodicBonds = []
    for bond in morphology['bond']:
        posn1 = np.array(morphology['position'][bond[1]]) + (np.array(morphology['image'][bond[1]]) *
                                        np.array([morphology['lx'], morphology['ly'], morphology['lz']]))
        posn2 = np.array(morphology['position'][bond[2]]) + (np.array(morphology['image'][bond[2]]) *
                                        np.array([morphology['lx'], morphology['ly'], morphology['lz']]))
        separation = helperFunctions.calculateSeparation(posn1, posn2)
        if separation >= morphology['lx'] / 2.0:
            print("Periodic bond found:", bond, "because separation =", separation, ">=", morphology['lx'] / 2.0)
            morphology = move_bonded_atoms(bond[1], morphology, bond_dict)
    return morphology


def zero_out_images(morphology):
    for atomID, image in enumerate(morphology['image']):
        if image != [0, 0, 0]:
            morphology['image'][atomID] = [0, 0, 0]
    return morphology


def getbond_dict(morphology):
    bond_dict = {atomID: [] for atomID, atomType in enumerate(morphology['type'])}
    for bond in morphology['bond']:
        #if bond[1] < bond[2]:
        bond_dict[bond[1]].append(bond[2])
        #else:
        bond_dict[bond[2]].append(bond[1])
    return bond_dict


def move_bonded_atoms(central_atom, morphology, bond_dict):
    for bonded_atom in bond_dict[central_atom]:
        atom1_posn = morphology['position'][central_atom]
        atom2_posn = morphology['position'][bonded_atom]
        #print("atom1:", central_atom, "posn =", atom1_posn, "; atom2:", bonded_atom, "posn =", atom2_posn)
        sep_vec = np.array(atom1_posn) - np.array(atom2_posn)
        moved = False
        for axis, value in enumerate(sep_vec):
            if value > morphology['lx'] / 2.0:
                morphology['position'][bonded_atom][axis] += morphology['lx']
                moved = True
            if value < -morphology['lx'] / 2.0:
                morphology['position'][bonded_atom][axis] -= morphology['lx']
                moved = True
        if moved:
            #print("Moved", bonded_atom, "to same box as", central_atom)
            #print("New Positions: atom1 =", morphology['position'][central_atom], 
            #       "atom2 =", morphology['position'][bonded_atom])
            morphology = move_bonded_atoms(bonded_atom, morphology, bond_dict)
        else:
            #print("Move was unnecessary")
            pass
    return morphology


def execute():
    list_of_files = sys.argv[1:]
    if len(list_of_files) < 1:
        print("No files requested to convert!")
        exit()
    for file_name in list_of_files:
        print("Fixing the images for", file_name + "...")
        morphology = helperFunctions.loadMorphologyXML(file_name)
        morphology = zero_out_images(morphology)
        bond_dict = getbond_dict(morphology)
        morphology = check_bonds(morphology, bond_dict)
        split_file_name = file_name.split('/')
        file_directory = '/'.join(split_file_name[:-1])
        if len(split_file_name) > 1:
            file_directory += '/'
        helperFunctions.writeMorphologyXML(morphology, file_directory + "imageFix_" + split_file_name[-1])


if __name__ == "__main__":
    execute()
