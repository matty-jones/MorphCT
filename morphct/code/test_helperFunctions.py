import unittest
import os
import sys
import copy
import numpy as np
import helperFunctions_modified as hf

test_morphology_dict = {'natoms': 4,
'dimensions': 3,
'lz': 40.0,
'dihedral': [['C10-C9-C2-C1', 0, 1, 2, 3]],
'yz': 0.0,
'improper': [['C10-C9-C2-C1', 0, 1, 2, 3]],
'xz': 0.0,
'charge': [0.0, 0.0, 0.0, 0.0],
'bond': [['C2-C9', 2, 1], ['C1-C2', 3, 2], ['C10-C9', 0, 1]],
'mass': [1.0, 1.0, 1.0, 1.0],
'diameter': [1.0, 1.0, 1.0, 1.0],
'time_step': 0,
'angle': [['C10-C9-C2', 0, 1, 2], ['C9-C2-C1', 1, 2, 3]],
'xy': 0.0,
'lx': 40.0,
'ly': 40.0,
'image': [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
'body': [-1, -1, -1, -1],
'position': [[-2.77009144, -4.36452856, 0.28093504], [-2.07113944, -3.15989656, 0.19025104],
             [-0.67753344, -3.26944556, 0.49046004], [-0.27401744, -4.56713156, 0.81986104]],
'type': ['C10', 'C9', 'C2', 'C1']}

test_FF_dict = {'lj': [['C1', 0.07, 3.55], ['C2', 0.07, 3.55], ['C9', 0.07, 3.55], ['C10', 0.07, 3.55]],
'improper': [['C10-C9-C2-C1', 126.32, -109.81, -19.738, -25.303, 28.53]],
'dihedral': [['C10-C9-C2-C1', 126.32, -109.81, -19.738, -25.303, 28.53]],
'angle': [['C9-C2-C1', 39.582, 1.97784], ['C10-C9-C2', 39.582, 1.97784],
          ['C1-C10-C9', 54.694, 2.27137], ['C2-C1-C10', 54.694, 2.27137]],
'bond': [['C1-C2', 1028.54, 1.37368], ['C10-C9', 1028.54, 1.37368],
         ['C2-C9', 906.2, 1.43277], ['C1-C10', 784.58, 1.45]],
'dpd': [['C1', 10.0, 2.0], ['C2', 10.0, 2.0], ['C9', 10.0, 2.0], ['C10', 10.0, 2.0]]}

class TestCommand(unittest.TestCase):
    file_created = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def fn_response(self, function, posn_args, kw_args):
        cmd = getattr(hf, function)
        if (posn_args is not None) and (kw_args is not None):
            if isinstance(posn_args, dict):
                return cmd(posn_args, **kw_args)
            else:
                return cmd(*posn_args, **kw_args)
        elif (posn_args is None):
            return cmd(**kw_args)
        elif (kw_args is None):
            if isinstance(posn_args, dict):
                return cmd(posn_args)
            else:
                return cmd(*posn_args)

    def compare_equal(self, function, expected, posn_args=None, kw_args=None):
        response = self.fn_response(function, posn_args, kw_args)
        #print("Expected", repr(expected), "got", repr(response))
        # Dictionary
        if isinstance(expected, (dict)):
            print(type(expected))
            for key, val in expected.items():
                print(type(val))
                if isinstance(val, (list, np.ndarray)):
                    # Array in dict (a la morphology_dictionary['position'])
                    for index, value in enumerate(val):
                        print(type(value))
                        if isinstance(value, (float)):
                            # Check that the answer is within 1E-4% of expected
                            difference = np.abs(value - response[key][index])
                            self.assertTrue(difference <= np.abs(1E-6 * value), msg="Expected " + repr(value)
                                            + " for key " + repr([key, index]) + ", but got "
                                            + repr(response[key][index]) + ", which is more than 1E-4% ("
                                            + repr(difference) + ") from expected.")
                        else:
                            self.assertEqual(response[key][index], value, msg="Expected " + repr(value) + " for key "
                                             + repr([key, index]) + ", but got " + repr(response[key][index])
                                             + " instead.")
                elif isinstance(val, (float)):
                    # Check that the answer is within 1E-4% of expected
                    difference = np.abs(val - response[key])
                    self.assertTrue(difference <= np.abs(1E-6 * val), msg="Expected " + repr(val) + " for key "
                                    + repr(key) + ", but got " + repr(response[key]) + ", which is more than 1E-4% ("
                                    + repr(difference) + ") from expected.")
                else:
                    self.assertEqual(response[key], val, msg="Expected " + repr(val) + " for key " + repr(key)
                                     + ", but got " + repr(response[key]) + " instead.")
        # Tensor
        elif isinstance(expected, (np.matrixlib.defmatrix.matrix)):
            for rowID in range(expected.shape[0]):
                for colID in range(expected.shape[1]):
                    val = float(expected[rowID, colID])
                    # Check that the answer is within 1E-4% of expected
                    difference = np.abs(val - response[rowID, colID])
                    self.assertTrue(difference <= np.abs(1E-6 * val), msg="Expected " + repr(val) + " for index "
                                    + repr([rowID, colID]) + ", but got " + repr(response[rowID, colID])
                                    + ", which is more than 1E-4% (" + repr(difference) + ") from expected.")
        # Vector
        elif isinstance(expected, (list, np.ndarray)):
            try:
                for index, val in enumerate(expected):
                    if isinstance(val, (float)):
                        # Check that the answer is within 1E-4% of expected
                        difference = np.abs(val - response[index])
                        self.assertTrue(difference <= np.abs(1E-6 * val), msg="Expected " + repr(val) + " for index "
                                        + repr(index) + ", but got " + repr(response[index])
                                        + ", which is more than 1E-4% (" + repr(difference) + ") from expected.")
                    else:
                        self.assertEqual(response[index], val, msg="Expected " + repr(val) + " for index "
                                         + repr(index) + ", but got " + repr(response[index]) + " instead.")
            except ValueError:
                # Actually a list of vectors
                for array_index, array in enumerate(expected):
                    for index, val in enumerate(array):
                        if isinstance(val, (float)):
                            # Check that the answer is within 1E-4% of expected
                            difference = np.abs(val - response[array_index][index])
                            self.assertTrue(difference <= np.abs(1E-6 * val), msg="Expected " + repr(val)
                                            + " for index " + repr([array_index, index]) + ", but got "
                                            + repr(response[array_index][index]) + ", which is more than 1E-4% ("
                                            + repr(difference) + ") from expected.")
                        else:
                            self.assertEqual(response[array_index][index], val, msg="Expected " + repr(val)
                                             + " for index " + repr([array_index, index]) + ", but got "
                                             + repr(response[array_index][index]) + " instead.")
        # Scalar
        else:
            if isinstance(response, (float)):
                # Check that the answer is within 1E-4% of expected
                difference = np.abs(expected - response)
                self.assertTrue(difference <= (1E-6 * expected), msg="Expected " + repr(response) + " but got "
                                + repr(response) + ", which is more than 1E-4% (" + repr(difference) + ")"
                                " from expected.")
            else:
                self.assertEqual(response, expected, msg="Expected " + repr(response) + " but got "
                                 + repr(response) + " instead.")

    def confirm_file_exists(self, function, expected, posn_args=None, kw_args=None):
        response = self.fn_response(function, posn_args, kw_args)
        directory = '/'.join(expected.split('/')[:-1])
        file_name = expected.split('/')[-1]
        files = os.listdir(directory)
        self.assertTrue(file_name in files)
        self.file_created = expected

    def setUp(self):
        #pass
        sys.stdout = None

    def tearDown(self):
        if self.file_created is not None:
            os.remove('./' + self.file_created)
        sys.stdout = sys.__stdout__


class TestMathematicalOperations(TestCommand):
#    def test_find_index(self):
#        function = "find_index"
#        # Exists once
#        self.compare_equal(function, [0], posn_args=["wobbey", "w"])
#        # Exists twice`
#        self.compare_equal(function, [2, 3], posn_args=["wobbey", "b"])
#        # Does not exist
#        self.compare_equal(function, None, posn_args=["wobbey", "z"])
#
#    def test_calculate_separation(self):
#        function = "calculate_separation"
#        # Check with integers
#        self.compare_equal(function, 1.0, posn_args=[[0, 0, 0], [1, 0, 0]])
#        # Check with floats
#        self.compare_equal(function, 8.395236745, posn_args=[[1.0, 2.0, 3.0], [4.0, 6.8, 9.2]])
#        # Check with numpy arrays
#        self.compare_equal(function, 23.17593795, posn_args=[np.array([-4.7, -2.9, -7.2]),
#                                                             np.array([6.8, 15.2, 1.59])]
#                          )
#
#    def test_calc_com(self):
#        function = "calc_com"
#        # Check COM with atom types
#        self.compare_equal(function, np.array([0.65551656, -0.26839579, 1.19994922]),
#                           posn_args=[[[9.785, 1.847, -0.12], [1.391, -8.481, 6.544], [-2.97, 5.524, -9.939],
#                                       [-3.374, -0.493, 9.536], [-1.09, -7.816, -1.172], [7.859, 5.693, 7.094],
#                                       [-2.048, -2.932, 6.502], [-5.084, -1.632, -5.404], [0.686, 7.565, -4.925],
#                                       [-1.699, -8.273, 5.301]]],
#                           kw_args={'list_of_atom_types': ['S', 'H', 'O', 'SI', 'S', 'O', 'SI', 'C', 'Si', 'O']}
#                              )
#        # Check COM with atom masses
#        self.compare_equal(function, np.array([-2.72260512,  4.30470998,  0.42963083]),
#                           posn_args = [[[8.849, 0.054, -8.658], [-0.988, 1.209, 1.179], [-9.316, 7.25, 6.749],
#                                         [-9.209, 5.239, -0.469], [8.19, 6.992, 1.705], [0.395, -0.266, -1.127],
#                                         [8.025, 9.111, 3.071], [-8.334, -0.709, -9.92], [2.898, -3.857, -3.284],
#                                         [-6.342, 5.16, -9.741]]],
#                           kw_args = {'list_of_atom_masses': [2.10132191, 46.51193374, 75.50176202, 66.9189811,
#                                                              12.89007787, 0.73334827, 57.03730002, 17.13006727,
#                                                              45.66377561, 24.76080088]}
#                              )
#
#    def test_find_axis(self):
#        function = "find_axis"
#        # Check axis with normalisation
#        self.compare_equal(function, np.array([1.08101166, -0.04656524, -0.03444642]),
#                           posn_args=[[-8.447, 2.508, -9.044], [7.966, 1.801, -9.567]]
#                          )
#        # Check axis without normalisation
#        self.compare_equal(function, np.array([-2.29, -6.992, 6.782]),
#                           posn_args=[[8.743, 7.955, -1.732], [6.453, 0.963, 5.05]],
#                           kw_args = {'normalise': False}
#                          )
#
#    def test_get_rotation_matrix(self):
#        function = "get_rotation_matrix"
#        self.compare_equal(function, np.matrix([[ 40.56407824, -70.10992864,  86.60426525],
#                                                [ 94.64591336,  31.92827421,  31.85780149],
#                                                [-58.80350275, -71.11516251,  27.01137555]]),
#                           posn_args=[[-8.34, -5.539, -0.324],
#                                      [0.219, -9.732, 8.726]]
#                          )
#
#    def test_parallel_sort(self):
#        function = "parallel_sort"
#        self.compare_equal(function, ((1, 2, 3, 4, 5, 7), ('one', 'two', 'three', 'four', 'five', 'seven')),
#                           posn_args=[[5, 4, 3, 7, 1, 2], ['five', 'four', 'three', 'seven', 'one', 'two']]
#                          )
#
#    def test_terminating_positions(self):
#        function = "get_terminating_positions"
#        # Add one hydrogen to perylene
#        current_atom = [2.626117, -2.652472, 0.304705]
#        bonded_neighbours = [[1.474839, -3.416274, 0.414203], [2.546457, -1.273502, 0.131632]]
#        output = [np.array([ 3.57512402, -3.12219452,  0.35314501])]
#        self.compare_equal(function, output, posn_args=[current_atom, bonded_neighbours, 1])
#        # Add two hydrogens to middle of alkyl sidechain
#        current_atom = [0.774749, 0.328046, 0.092399]
#        bonded_neighbours = [[-0.294037, -0.771723, 0.088898], [0.169728, 1.684461, -0.288625]]
#        output = [np.array([ 1.19729525,  0.39736037,  1.06206399]),
#                  np.array([ 1.53157525,  0.07780889, -0.60631082])]
#        self.compare_equal(function, output, posn_args=[current_atom, bonded_neighbours, 2])
#        # Add three hydrogens to end of alkyl sidechain
#        current_atom = [0.640141, 4.135119, -0.663760]
#        bonded_neighbours = [[1.238487, 2.783469, -0.284985]]
#        output = [np.array([-0.32155015,  4.23087175, -0.22835037]),
#                  np.array([ 1.26746188,  4.91197293, -0.30799815]),
#                  np.array([ 0.55827332,  4.20278437, -1.71842539])]
#        self.compare_equal(function, output, posn_args=[current_atom, bonded_neighbours, 3])
#
#
#class TestFileManipHelperFunctions(TestCommand):
#    def test_write_CSV(self):
#        function = "write_CSV"
#        file_name = "./testAssets/test.csv"
#        self.confirm_file_exists(function, file_name,
#                                 posn_args=[file_name, [['el1_1', 'el1_2'], ['row2_1', 'row2_2']]],
#                                )
#
#    def test_load_morphology_xml(self):
#        function = "load_morphology_XML"
#        input_xml = "./testAssets/test_input.xml"
#        output_dictionary = copy.deepcopy(test_morphology_dict)
#        self.compare_equal(function, output_dictionary, posn_args=[input_xml])
#
#    def test_load_FF_xml(self):
#        function = "load_FF_XML"
#        input_xml = "./testAssets/test_FF.xml"
#        output_dictionary = copy.deepcopy(test_FF_dict)
#        self.compare_equal(function, output_dictionary, posn_args=[input_xml])
#
#    def test_write_morphology_xml(self):
#        function = "write_morphology_XML"
#        input_dictionary = copy.deepcopy(test_morphology_dict)
#        output_xml = "./testAssets/test_output.xml"
#        self.confirm_file_exists(function, output_xml, posn_args=[input_dictionary, output_xml])
#
#    def test_write_XYZ_file(self):
#        function = "write_XYZ_file"
#        input_dictionary = copy.deepcopy(test_morphology_dict)
#        output_xml = "./testAssets/test_output.xyz"
#        self.confirm_file_exists(function, output_xml, posn_args=[input_dictionary, output_xml])
#
#
#class TestDictManipHelperFunctions(TestCommand):
#    def test_add_unwrapped_positions(self):
#        function = "add_unwrapped_positions"
#        input_dictionary = {'lz': 5.0, 'ly': 5.0, 'image': [[9, -8, 7], [-6, 5, -4], [3, -2, 1]],
#                            'position': [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]], 'lx': 5.0}
#        output_dictionary = {'unwrapped_position': [[46.0, -38.0, 38.0], [-26.0, 30.0, -14.0], [22.0, -2.0, 14.0]],
#                             'lz': 5.0, 'ly': 5.0, 'image': [[9, -8, 7], [-6, 5, -4], [3, -2, 1]],
#                             'position': [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]], 'lx': 5.0}
#        self.compare_equal(function, output_dictionary, posn_args=input_dictionary)
#
#    def test_replace_wrapped_positions(self):
#        function = "replace_wrapped_positions"
#        input_dictionary = {'unwrapped_position': [[46.0, -38.0, 38.0], [-26.0, 30.0, -14.0], [22.0, -2.0, 14.0]],
#                            'lz': 5.0, 'ly': 5.0, 'image': [[9, -8, 7], [-6, 5, -4], [3, -2, 1]],
#                            'position': [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]], 'lx': 5.0}
#        output_dictionary = {'unwrapped_position': [[46.0, -38.0, 38.0], [-26.0, 30.0, -14.0], [22.0, -2.0, 14.0]],
#                             'lz': 5.0, 'ly': 5.0, 'image': [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
#                             'position': [[46.0, -38.0, 38.0], [-26.0, 30.0, -14.0], [22.0, -2.0, 14.0]], 'lx': 5.0}
#        self.compare_equal(function, output_dictionary, posn_args=input_dictionary)
#
#    def test_add_masses(self):
#        function = "add_masses"
#        input_dictionary = {'type': ['H', 'C', 'O', 'H', 'H', 'S', 'N', 'C', 'H', 'N']}
#        output_dictionary = {'mass': [1.007825, 12.0, 15.994914, 1.007825, 1.007825, 31.972071, 14.003074, 12.0,
#                                      1.007825, 14.003074],
#                             'type': ['H', 'C', 'O', 'H', 'H', 'S', 'N', 'C', 'H', 'N']}
#        self.compare_equal(function, output_dictionary, posn_args=input_dictionary)
#
#    def test_add_diameters(self):
#        function = "add_diameters"
#        input_dictionary = {'type': ['H', 'C', 'O', 'H', 'H', 'S', 'N', 'C', 'H', 'N']}
#        output_dictionary = {'diameter': [1.2, 1.7, 1.52, 1.2, 1.2, 1.8, 1.55, 1.7, 1.2, 1.55],
#                             'type': ['H', 'C', 'O', 'H', 'H', 'S', 'N', 'C', 'H', 'N']}
#        self.compare_equal(function, output_dictionary, posn_args=input_dictionary)
#
#    def test_check_constraint_names(self):
#        function = "check_constraint_names"
#        # First copy the original dictionary
#        input_dictionary = copy.deepcopy(test_morphology_dict)
#        # Now change some of the constraints
#        input_dictionary['bond'][0][0] = 'test1'
#        input_dictionary['angle'][0][0] = 'test2'
#        input_dictionary['dihedral'][0][0] = 'test3'
#        input_dictionary['improper'][0][0] = 'test4'
#        output_dictionary = copy.deepcopy(test_morphology_dict)
#        self.compare_equal(function, output_dictionary, posn_args=[input_dictionary])
#
#    def test_increment_atom_IDs(self):
#        function = "increment_atom_IDs"
#        input_dictionary = copy.deepcopy(test_morphology_dict)
#        # Create the ghost_dictionary
#        ghost_input = copy.deepcopy(test_morphology_dict)
#        ghost_input['bond'][0][1] = '_2'
#        ghost_input['bond'][1][2] = '_2'
#        ghost_input['bond'][2][1] = '_0'
#        ghost_input['bond'][2][2] = '_1'
#        # Create the output_dictionary
#        output_dictionary = copy.deepcopy(input_dictionary)
#        output_dictionary['bond'] = [['C2-C9', 27, 26], ['C1-C2', 28, 27], ['C10-C9', 25, 26]]
#        output_dictionary['angle'] = [['C10-C9-C2', 25, 26, 27], ['C9-C2-C1', 26, 27, 28]]
#        output_dictionary['dihedral'] = [['C10-C9-C2-C1', 25, 26, 27, 28]]
#        output_dictionary['improper'] = [['C10-C9-C2-C1', 25, 26, 27, 28]]
#        # Increment without modifying the ghost dictionary
#        self.compare_equal(function, [output_dictionary], posn_args=[input_dictionary, ghost_input, 25])
#        # Create the ghost output dictionary
#        ghost_output = copy.deepcopy(test_morphology_dict)
#        ghost_output['bond'][0][1] = 27
#        ghost_output['bond'][1][2] = 27
#        ghost_output['bond'][2][1] = 25
#        ghost_output['bond'][2][2] = 26
#        # Increment and modify the ghost dictionary
#        self.compare_equal(function, [output_dictionary, ghost_output], posn_args=[input_dictionary,
#                                                                                   ghost_input, 25],
#                           kw_args={'modify_ghost_dictionary': True})
#
#    def test_scale(self):
#        function = "scale"
#        input_dictionary = copy.deepcopy(test_morphology_dict)
#        output_dictionary = copy.deepcopy(test_morphology_dict)
#        output_dictionary['position'] = [[-5.5401828799999997, -8.7290571200000002, 0.56187008000000005],
#                                         [-4.1422788800000001, -6.3197931199999999, 0.38050208000000002],
#                                         [-1.3550668800000001, -6.5388911199999997, 0.98092007999999997],
#                                         [-0.54803488, -9.13426312, 1.6397220800000001]]
#        output_dictionary['lx'] = 80.0
#        output_dictionary['ly'] = 80.0
#        output_dictionary['lz'] = 80.0
#        self.compare_equal(function, output_dictionary, posn_args=[input_dictionary, 2])

    def test_rotate(self):
        function = "rotate"
        # Rotate 37 degrees around origin
        theta = 0.6457718232379019  # 37 degrees
        input_dictionary = copy.deepcopy(test_morphology_dict)
        output_dictionary = copy.deepcopy(test_morphology_dict)
        output_dictionary['position'] = [[0.41434546632213598, -5.1527501367284678, 0.28093504000000002],
                                         [0.24758771837101312, -3.7700484309270061, 0.19025104000000004],
                                         [1.4264991909572644, -3.0188451252623438, 0.49046003999999999],
                                         [2.5297283275635971, -3.8123812548713532, 0.81986104000000004]]
        output_dictionary['unwrapped_position'] = output_dictionary['position']
        self.compare_equal(function, output_dictionary, posn_args=[input_dictionary, theta])
        # Rotate 37 degrees around [1, 1, 2]
        input_dictionary = copy.deepcopy(test_morphology_dict)
        output_dictionary = copy.deepcopy(test_morphology_dict)
        output_dictionary['position'] = [[1.2175249794268916, -5.5532006699278078, 0.28093504000000002],
                                         [1.0507672314757688, -4.170498964126347, 0.19025104000000004],
                                         [2.2296787040620196, -3.4192956584616852, 0.49046003999999999],
                                         [3.3329078406683523, -4.2128317880706945, 0.81986104000000004]]
        output_dictionary['unwrapped_position'] = output_dictionary['position']
        self.compare_equal(function, output_dictionary, posn_args=[input_dictionary, theta],
                           kw_args={'rotate_around_point': [1, 1, 2]})
        # Rotate 37 degrees around [1, 1, 2] and axis [-1, -5, -3] (normalized)
        input_dictionary = copy.deepcopy(test_morphology_dict)
        output_dictionary = copy.deepcopy(test_morphology_dict)
        output_dictionary['position'] = [[-2.9793708141777939, -3.3370294782402148, -1.3618036382070449],
                                         [-1.9703050205356663, -2.4119510409454508, -1.0899362982456937],
                                         [-1.033399740537023, -2.8439568033438838, -0.10006578758118634],
                                         [-1.3040248954050759, -4.1165836731432366, 0.41228371370708689]]
        output_dictionary['unwrapped_position'] = output_dictionary['position']
        self.compare_equal(function, output_dictionary, posn_args=[input_dictionary, theta],
                           kw_args={'rotate_around_point': [1, 1, 2],
                                    'rotate_around_axis': [-0.16903085, -0.84515425, -0.50709255]})
        # Rotate 37 degrees around [1, 1, 2] and axis [-1, -5, -3] (not normalized)
        input_dictionary = copy.deepcopy(test_morphology_dict)
        self.compare_equal(function, output_dictionary, posn_args=[input_dictionary, theta],
                           kw_args={'rotate_around_point': [1, 1, 2],
                                    'rotate_around_axis': [-1, -5, -3]})



if __name__ == "__main__":
    unittest.main()
