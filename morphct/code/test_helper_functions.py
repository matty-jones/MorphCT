import unittest
import os
import numpy as np
import helperFunctions_modified as hf

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
            for key, val in expected.items():
                if isinstance(val, (float)):
                    # Check that the answer is within 1E-4% of expected
                    difference = np.abs(val - response[key])
                    self.assertTrue(difference < np.abs(1E-6 * val))
                else:
                    self.assertEqual(response[key], val)
        # Tensor
        elif isinstance(expected, (np.matrixlib.defmatrix.matrix)):
            for rowID in range(expected.shape[0]):
                for colID in range(expected.shape[1]):
                    val = float(expected[rowID, colID])
                    # Check that the answer is within 1E-4% of expected
                    difference = np.abs(val - response[rowID, colID])
                    self.assertTrue(difference < np.abs(1E-6 * val))
        # Vector
        elif isinstance(expected, (list, np.ndarray)):
            for index, val in enumerate(expected):
                if isinstance(val, (float)):
                    # Check that the answer is within 1E-4% of expected
                    difference = np.abs(val - response[index])
                    self.assertTrue(difference < np.abs(1E-6 * val))
                else:
                    self.assertEqual(response[index], val)
        # Scalar
        else:
            if isinstance(response, (float)):
                # Check that the answer is within 1E-4% of expected
                difference = np.abs(expected - response)
                self.assertTrue(difference < (1E-6 * expected))
            else:
                self.assertEqual(response, expected)

    def confirm_file_exists(self, function, expected, posn_args=None, kw_args=None):
        response = self.fn_response(function, posn_args, kw_args)
        files = os.listdir('./')
        self.assertTrue(expected in files)
        self.file_created = expected

    def setUp(self):
        pass

    def tearDown(self):
        if self.file_created is not None:
            os.remove('./' + self.file_created)


class TestBasicHelperFunctions(TestCommand):
    def test_find_index(self):
        function = "find_index"
        # Exists once
        self.compare_equal(function, [0], posn_args=["wobbey", "w"])
        # Exists twice`
        self.compare_equal(function, [2, 3], posn_args=["wobbey", "b"])
        # Does not exist
        self.compare_equal(function, None, posn_args=["wobbey", "z"])

    def test_calculate_separation(self):
        function = "calculate_separation"
        # Check with integers
        self.compare_equal(function, 1.0, posn_args=[[0, 0, 0], [1, 0, 0]])
        # Check with floats
        self.compare_equal(function, 8.395236745, posn_args=[[1.0, 2.0, 3.0], [4.0, 6.8, 9.2]])
        # Check with numpy arrays
        self.compare_equal(function, 23.17593795, posn_args=[np.array([-4.7, -2.9, -7.2]),
                                                             np.array([6.8, 15.2, 1.59])]
                          )

    def test_calc_com(self):
        function = "calc_com"
        # Check COM with atom types
        self.compare_equal(function, np.array([0.65551656, -0.26839579, 1.19994922]),
                           posn_args=[[[9.785, 1.847, -0.12], [1.391, -8.481, 6.544], [-2.97, 5.524, -9.939],
                                       [-3.374, -0.493, 9.536], [-1.09, -7.816, -1.172], [7.859, 5.693, 7.094],
                                       [-2.048, -2.932, 6.502], [-5.084, -1.632, -5.404], [0.686, 7.565, -4.925],
                                       [-1.699, -8.273, 5.301]]],
                           kw_args={'list_of_atom_types': ['S', 'H', 'O', 'SI', 'S', 'O', 'SI', 'C', 'Si', 'O']}
                              )
        # Check COM with atom masses
        self.compare_equal(function, np.array([-2.72260512,  4.30470998,  0.42963083]),
                           posn_args = [[[8.849, 0.054, -8.658], [-0.988, 1.209, 1.179], [-9.316, 7.25, 6.749],
                                         [-9.209, 5.239, -0.469], [8.19, 6.992, 1.705], [0.395, -0.266, -1.127],
                                         [8.025, 9.111, 3.071], [-8.334, -0.709, -9.92], [2.898, -3.857, -3.284],
                                         [-6.342, 5.16, -9.741]]],
                           kw_args = {'list_of_atom_masses': [2.10132191, 46.51193374, 75.50176202, 66.9189811,
                                                              12.89007787, 0.73334827, 57.03730002, 17.13006727,
                                                              45.66377561, 24.76080088]}
                              )

    def test_find_axis(self):
        function = "find_axis"
        # Check axis with normalisation
        self.compare_equal(function, np.array([1.08101166, -0.04656524, -0.03444642]),
                           posn_args=[[-8.447, 2.508, -9.044], [7.966, 1.801, -9.567]]
                          )
        # Check axis without normalisation
        self.compare_equal(function, np.array([-2.29, -6.992, 6.782]),
                           posn_args=[[8.743, 7.955, -1.732], [6.453, 0.963, 5.05]],
                           kw_args = {'normalise': False}
                          )

    def test_get_rotation_matrix(self):
        function = "get_rotation_matrix"
        self.compare_equal(function, np.matrix([[ 40.56407824, -70.10992864,  86.60426525],
                                                [ 94.64591336,  31.92827421,  31.85780149],
                                                [-58.80350275, -71.11516251,  27.01137555]]),
                           posn_args=[[-8.34, -5.539, -0.324],
                                      [0.219, -9.732, 8.726]]
                          )

    def test_parallel_sort(self):
        function = "parallel_sort"
        self.compare_equal(function, ((1, 2, 3, 4, 5, 7), ('one', 'two', 'three', 'four', 'five', 'seven')),
                           posn_args=[[5, 4, 3, 7, 1, 2], ['five', 'four', 'three', 'seven', 'one', 'two']]
                          )

    def test_write_CSV(self):
        function = "write_CSV"
        file_name = "test.csv"
        self.confirm_file_exists(function, file_name,
                                 posn_args=[file_name, [['el1_1', 'el1_2'], ['row2_1', 'row2_2']]],
                                 kw_args = {'test': True}
                                )


class TestDictManipHelperFunctions(TestCommand):
    def test_add_unwrapped_positions(self):
        function = "add_unwrapped_positions"
        input_dictionary = {'lz': 5.0, 'ly': 5.0, 'image': [[9, -8, 7], [-6, 5, -4], [3, -2, 1]],
                            'position': [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]], 'lx': 5.0}
        output_dictionary = {'unwrapped_position': [[46.0, -38.0, 38.0], [-26.0, 30.0, -14.0], [22.0, -2.0, 14.0]],
                             'lz': 5.0, 'ly': 5.0, 'image': [[9, -8, 7], [-6, 5, -4], [3, -2, 1]],
                             'position': [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]], 'lx': 5.0}
        self.compare_equal(function, output_dictionary, posn_args=input_dictionary)

    def test_replace_wrapped_positions(self):
        function = "replace_wrapped_positions"
        input_dictionary = {'unwrapped_position': [[46.0, -38.0, 38.0], [-26.0, 30.0, -14.0], [22.0, -2.0, 14.0]],
                            'lz': 5.0, 'ly': 5.0, 'image': [[9, -8, 7], [-6, 5, -4], [3, -2, 1]],
                            'position': [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]], 'lx': 5.0}
        output_dictionary = {'unwrapped_position': [[46.0, -38.0, 38.0], [-26.0, 30.0, -14.0], [22.0, -2.0, 14.0]],
                             'lz': 5.0, 'ly': 5.0, 'image': [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
                             'position': [[46.0, -38.0, 38.0], [-26.0, 30.0, -14.0], [22.0, -2.0, 14.0]], 'lx': 5.0}
        self.compare_equal(function, output_dictionary, posn_args=input_dictionary)

    def test_add_masses(self):
        function = "add_masses"
        input_dictionary = {'type': ['H', 'C', 'O', 'H', 'H', 'S', 'N', 'C', 'H', 'N']}
        output_dictionary = {'mass': [1.007825, 12.0, 15.994914, 1.007825, 1.007825, 31.972071, 14.003074, 12.0,
                                      1.007825, 14.003074], 'type': ['H', 'C', 'O', 'H', 'H', 'S', 'N', 'C', 'H', 'N']}
        self.compare_equal(function, output_dictionary, posn_args=input_dictionary)

    def test_add_diameters(self):
        function = "add_diameters"
        input_dictionary = {'type': ['H', 'C', 'O', 'H', 'H', 'S', 'N', 'C', 'H', 'N']}
        output_dictionary = {'diameter': [1.2, 1.7, 1.52, 1.2, 1.2, 1.8, 1.55, 1.7, 1.2, 1.55],
                             'type': ['H', 'C', 'O', 'H', 'H', 'S', 'N', 'C', 'H', 'N']}
        self.compare_equal(function, output_dictionary, posn_args=input_dictionary)


if __name__ == "__main__":
    unittest.main()
