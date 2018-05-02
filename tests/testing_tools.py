import os
import sys
import shutil
import numpy as np
import time as T
from morphct.code import helper_functions as hf
from morphct.definitions import TEST_ROOT


class TestCommand(object):
    file_created = None

    def fn_response(self, function, posn_args, kw_args):
        cmd = getattr(hf, function)
        if (posn_args is not None) and (kw_args is not None):
            if isinstance(posn_args, dict):
                return cmd(posn_args, **kw_args)
            else:
                return cmd(*posn_args, **kw_args)
        elif (posn_args is None) and (kw_args is not None):
            return cmd(**kw_args)
        elif (kw_args is None) and (posn_args is not None):
            if isinstance(posn_args, dict):
                return cmd(posn_args)
            else:
                return cmd(*posn_args)
        else:
            # No args
            return cmd()

    def compare_lt(self, object1, object2):
        assert object2 > object1, ("Expected object1 (" + repr(object1)
                                   + ") to be less than object2 (" + repr(object2) + ").")

    def compare_equal(self, expected, response=None, function=None, posn_args=None, kw_args=None):
        if function is not None:
            response = self.fn_response(function, posn_args, kw_args)
        # Dictionary
        if isinstance(expected, (dict)):
            for key, val in expected.items():
                if isinstance(val, (list, np.ndarray)):
                    # Array in dict (a la morphology_dictionary['position'])
                    for index, value in enumerate(val):
                        # Coordinates in array (a la morphology_dictionary['position'])
                        if isinstance(value, (list, np.ndarray)):
                            for index2, value2 in enumerate(value):
                                if isinstance(value2, (float)):
                                    # Check that the answer is within 1E-4% of expected
                                    difference = np.abs(value2 - response[key][index][index2])
                                    assert difference <= np.abs(1E-6 * value2),\
                                        ("Expected " + repr(value2)
                                         + " for key " + repr([key, index, index2]) + ", but got "
                                         + repr(response[key][index][index2]) + ", which is more than "
                                         + "1E-4% (" + repr(difference) + ") from expected.")
                                else:
                                    assert response[key][index][index2] == value2,\
                                        ("Expected "
                                         + repr(value2) + " for key " + repr([key, index]) + ", but got "
                                         + repr(response[key][index][index2]) + " instead.")
                        elif isinstance(value, (float)):
                            # Check that the answer is within 1E-4% of expected
                            difference = np.abs(value - response[key][index])
                            assert difference <= np.abs(1E-6 * value),\
                                ("Expected " + repr(value)
                                 + " for key " + repr([key, index]) + ", but got "
                                 + repr(response[key][index]) + ", which is more than 1E-4% ("
                                 + repr(difference) + ") from expected.")
                        else:
                            assert response[key][index] == value,\
                                ("Expected " + repr(value) + " for key "
                                 + repr([key, index]) + ", but got " + repr(response[key][index])
                                 + " instead.")
                elif isinstance(val, (float)):
                    # Check that the answer is within 1E-4% of expected
                    difference = np.abs(val - response[key])
                    assert difference <= np.abs(1E-6 * val),\
                        ("Expected " + repr(val) + " for key "
                         + repr(key) + ", but got " + repr(response[key]) + ", which is more than 1E-4% ("
                         + repr(difference) + ") from expected.")
                else:
                    assert response[key] == val, ("Expected " + repr(val) + " for key " + repr(key)
                                                  + ", but got " + repr(response[key]) + " instead.")
        # Tensor
        elif isinstance(expected, (np.matrixlib.defmatrix.matrix)):
            for rowID in range(expected.shape[0]):
                for colID in range(expected.shape[1]):
                    val = float(expected[rowID, colID])
                    # Check that the answer is within 1E-4% of expected
                    difference = np.abs(val - response[rowID, colID])
                    assert difference <= np.abs(1E-6 * val),\
                        ("Expected " + repr(val) + " for index "
                         + repr([rowID, colID]) + ", but got " + repr(response[rowID, colID])
                         + ", which is more than 1E-4% (" + repr(difference) + ") from expected.")
        # Vector
        elif isinstance(expected, (list, np.ndarray)):
            try:
                for index, val in enumerate(expected):
                    if isinstance(val, (float)):
                        # Check that the answer is within 1E-4% of expected
                        difference = np.abs(val - response[index])
                        assert difference <= np.abs(1E-6 * val),\
                            ("Expected " + repr(val) + " for index "
                             + repr(index) + ", but got " + repr(response[index])
                             + ", which is more than 1E-4% (" + repr(difference) + ") from expected.")
                    else:
                        assert response[index] == val,\
                            ("Expected " + repr(val) + " for index "
                             + repr(index) + ", but got " + repr(response[index]) + " instead.")
            except ValueError:
                # Actually a list of vectors
                for array_index, array in enumerate(expected):
                    for index, val in enumerate(array):
                        if isinstance(val, (float)):
                            # Check that the answer is within 1E-4% of expected
                            difference = np.abs(val - response[array_index][index])
                            assert difference <= np.abs(1E-6 * val),\
                                ("Expected " + repr(val)
                                 + " for index " + repr([array_index, index]) + ", but got "
                                 + repr(response[array_index][index]) + ", which is more than 1E-4% ("
                                 + repr(difference) + ") from expected.")
                        else:
                            assert response[array_index][index] == val,\
                                ("Expected " + repr(val)
                                 + " for index " + repr([array_index, index]) + ", but got "
                                 + repr(response[array_index][index]) + " instead.")
        # Scalar
        else:
            if isinstance(response, (float)):
                # Check that the answer is within 1E-4% of expected
                difference = np.abs(expected - response)
                assert difference <= (1E-6 * expected),\
                    ("Expected " + repr(expected) + " but got "
                     + repr(response) + ", which is more than 1E-4% (" + repr(difference) + ")"
                     " from expected.")
            else:
                assert response == expected, ("Expected " + repr(expected) + " but got "
                                              + repr(response) + " instead.")

    def confirm_file_exists(self, expected, response=None, function=None, posn_args=None, kw_args=None):
        if function is not None:
            _ = self.fn_response(function, posn_args, kw_args)
        directory = '/'.join(expected.split('/')[:-1])
        file_name = expected.split('/')[-1]
        files = os.listdir(directory)
        assert file_name in files, ("Expected the file " + str(file_name) + " to exist, but it doesn't.")

    def setup_method(self):
        np.random.seed(929292929)
        sys.stdout = None
