import csv
import os
import sys
import shutil
import numpy as np
import time as T
import scipy.sparse
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
        assert object2 > object1, "".join(
            [
                "Expected object1 (",
                repr(object1),
                ") to be less than object2 (",
                repr(object2),
                ").",
            ]
        )

    def compare_equal(
        self,
        expected,
        response=None,
        function=None,
        posn_args=None,
        kw_args=None,
        dict_key=None,
    ):
        if function is not None:
            response = self.fn_response(function, posn_args, kw_args)
        if isinstance(expected, (float, int)):
            self.check_scalar(response, expected, dict_key)
        elif isinstance(expected, (np.matrixlib.defmatrix.matrix)):
            self.check_matrix(response, expected, dict_key)
        elif isinstance(expected, (scipy.sparse.lil.lil_matrix)):
            self.check_matrix(response.toarray(), expected.toarray(), dict_key)
        elif isinstance(expected, (list, np.ndarray)):
            self.check_array(response, expected, dict_key)
        elif isinstance(expected, (dict)):
            self.check_dictionary(response, expected, dict_key)
        elif isinstance(expected, (oc.chromophore)):
            self.check_dictionary(response.__dict__, expected.__dict__, dict_key)
        else:
            self.check_identical(response, expected, dict_key)

    def check_scalar(self, response, expected, dict_key=None):
        if dict_key is not None:
            assert_string = "".join(
                [
                    "Expected ",
                    repr(expected),
                    " for key ",
                    repr(dict_key),
                    " but got ",
                    repr(response),
                    ", which is more than 1E-4% from expected.",
                ]
            )
        else:
            assert_string = "".join(
                [
                    "Expected ",
                    repr(expected),
                    " but got ",
                    repr(response),
                    ", which is more than 1E-4% from expected.",
                ]
            )
        # Check that the answer is within 1E-4% of expected
        assert np.isclose(response, expected, 1E-6), assert_string

    def check_identical(self, response, expected, dict_key=None):
        if dict_key is not None:
            assert_string = " ".join(
                [
                    "Expected ",
                    repr(expected),
                    " for key ",
                    repr(dict_key),
                    " but got ",
                    repr(response),
                    " instead.",
                ]
            )
        else:
            assert_string = " ".join(
                ["Expected", repr(expected), "but got", repr(response), "instead."]
            )
        assert response == expected, assert_string

    def check_array(self, response, expected, dict_key=None):
        for expected_index, expected_value in enumerate(expected):
            self.compare_equal(
                response[expected_index], expected_value, dict_key=dict_key
            )

    def check_matrix(self, response, expected, dict_key=None):
        for expected_row_id in range(expected.shape[0]):
            for expected_col_id in range(expected.shape[1]):
                self.compare_equal(
                    response[expected_row_id, expected_col_id],
                    expected[expected_row_id, expected_col_id],
                    dict_key=dict_key,
                )

    def check_dictionary(self, response, expected, dict_key=None):
        for expected_key, expected_val in expected.items():
            if dict_key is None:
                check_dict_key = repr(expected_key)
            else:
                check_dict_key = " | ".join([repr(dict_key), repr(expected_key)])
            self.compare_equal(
                response[expected_key], expected_val, dict_key=check_dict_key
            )

    def confirm_file_exists(
        self,
        expected,
        response=None,
        function=None,
        posn_args=None,
        kw_args=None,
        negate=False,
    ):
        if function is not None:
            _ = self.fn_response(function, posn_args, kw_args)
        (directory, file_name) = os.path.split(expected)
        files = os.listdir(directory)
        if negate is False:
            assert file_name in files, "".join(
                [
                    "Expected the file ",
                    str(file_name),
                    " to exist in ",
                    str(directory),
                    ", but it doesn't.",
                ]
            )
        else:
            assert file_name not in files, "".join(
                [
                    "Expected the file ",
                    str(file_name),
                    " to not exist in ",
                    str(directory),
                    ", but it does.",
                ]
            )
