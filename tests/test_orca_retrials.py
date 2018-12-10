import csv
import os
import pytest
import shutil
import subprocess
from morphct import run_MorphCT
from morphct.definitions import TEST_ROOT
from testing_tools import TestCommand
from morphct.code import helper_functions as hf
from morphct.code import transfer_integrals as ti


class dummy_chromophore:
    def __init__(self):
        self.AAIDs = [0, 1]


@pytest.fixture(scope="module", params=[1, 3, 6, 9, 12, 15, 18])
def run_simulation(request):
    fail_count = request.param

    # ---==============================================---
    # ---=============== Setup Prereqs ================---
    # ---==============================================---

    output_dir = os.path.join(TEST_ROOT, "output_OR")
    orca_inp_dir = os.path.join(output_dir, "chromophores", "input_orca", "single")
    orca_out_dir = os.path.join(output_dir, "chromophores", "output_orca", "single")
    try:
        shutil.rmtree(output_dir)
    except OSError:
        pass
    os.makedirs(orca_inp_dir)
    os.makedirs(orca_out_dir)
    asset_name = os.path.join(TEST_ROOT, "assets", "orca_retrials", "00000.inp")
    shutil.copy(asset_name, orca_inp_dir)
    file_name = os.path.join("single", os.path.split(asset_name)[1])

    # Create dummy inputs for rerun_fails
    # First the failed_chromo_files dictionary that keeps track of the fail count and
    # the chromo_ID
    failed_chromo_files = {file_name: [fail_count, 0]}
    # Then, the parameter dictionary (used to get the orca output directory and the
    # proc_IDs
    parameter_dict = {
        "proc_IDs": [0],
        "output_orca_directory": output_dir,
        "output_morphology_directory": output_dir,
    }
    # Finally a blank chromophore class (used only to output AAIDs after 18 fails)
    chromophore_list = [dummy_chromophore()]

    ti.rerun_fails(failed_chromo_files, parameter_dict, chromophore_list)
    modified_file = os.path.join(os.path.split(orca_inp_dir)[:-1][0], file_name)
    return [fail_count, modified_file]


# ---==============================================---
# ---================= Run Tests ==================---
# ---==============================================---


class TestCompareOutputs(TestCommand):
    def test_check_output_exists(self, run_simulation):
        # No output is made for 18, so just return True
        if run_simulation[0] == 18:
            return True
        self.confirm_file_exists(
            os.path.join(
                TEST_ROOT,
                "output_OR",
                "chromophores",
                "output_orca",
                "single",
                "00000.out",
            )
        )

    def test_check_orca_ran_correctly(self, run_simulation):
        # No output is made for 18, so just return True
        if run_simulation[0] == 18:
            return True
        output_file = os.path.join(
            TEST_ROOT, "output_OR", "chromophores", "output_orca", "single", "00000.out"
        )
        with open(output_file, "r") as output_fh:
            output_lines = output_fh.readlines()
        SCF_success = False
        terminated_normally = False
        for line in output_lines:
            if SCF_success and terminated_normally:
                return True
            if "SUCCESS" in line:
                SCF_success = True
            if "ORCA TERMINATED NORMALLY" in line:
                terminated_normally = True
        return False

    def test_check_input_was_modified(self, run_simulation):
        fail_count = run_simulation[0]
        modified_inp_name = run_simulation[1]
        check_inp_name = os.path.join(
            TEST_ROOT, "assets", "orca_retrials", "00000_{:02d}.inp".format(fail_count)
        )
        with open(check_inp_name, "r") as expected_file:
            expected_lines = expected_file.readlines()
        with open(modified_inp_name, "r") as results_file:
            results_lines = results_file.readlines()
        self.compare_equal(expected_lines, response=results_lines)


def teardown_module():
    shutil.rmtree(os.path.join(TEST_ROOT, "output_OR"))


if __name__ == "__main__":

    class parameters:
        def __init__(self, param):
            self.param = param

    run_simulation(parameters(1))
