import os
import sys
import shutil
import numpy as np
import time as T
from morphct.code import helper_functions as hf
from morphct.code import obtain_chromophores as oc
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
        if isinstance(expected, (float, int)):
            self.check_scalar(response, expected)
        elif isinstance(expected, (np.matrixlib.defmatrix.matrix)):
            self.check_matrix(response, expected)
        elif isinstance(expected, (list, np.ndarray)):
            self.check_array(response, expected)
        elif isinstance(expected, (dict)):
            self.check_dictionary(response, expected)
        elif isinstance(expected, (oc.chromophore)):
            self.check_dictionary(response.__dict__, expected.__dict__)
        else:
            self.check_identical(response, expected)

    def check_scalar(self, response, expected):
        # Check that the answer is within 1E-4% of expected
        difference = np.abs(expected - response)
        assert difference <= np.abs(1E-6 * expected),\
            ("Expected " + repr(expected) + " but got "
             + repr(response) + ", which is more than 1E-4% (" + repr(difference) + ")"
             " from expected.")

    def check_identical(self, response, expected):
        assert response == expected, ("Expected " + repr(expected) + " but got "
                                      + repr(response) + " instead.")

    def check_array(self, response, expected):
        for expected_index, expected_value in enumerate(expected):
            self.compare_equal(response[expected_index], expected_value)

    def check_matrix(self, response, expected):
        for expected_row_id in range(expected.shape[0]):
            for expected_col_id in range(expected.shape[1]):
                self.compare_equal(response[expected_row_id, expected_col_id],
                                    expected[expected_row_id, expected_col_id])

    def check_dictionary(self, response, expected):
        for expected_key, expected_val in expected.items():
            self.compare_equal(response[expected_key], expected_val)

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
