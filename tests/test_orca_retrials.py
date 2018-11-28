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


@pytest.fixture(scope="module", params = [1, 3, 6, 9, 12, 15, 18])
def run_simulation(request):
    fail_count = request.param

    # ---==============================================---
    # ---=============== Setup Prereqs ================---
    # ---==============================================---

    output_dir = os.path.join(TEST_ROOT, "output_OR")
    orca_inp_dir = os.path.join(output_dir, "chromophores", "input_orca", "single")
    orca_out_dir = os.path.join(output_dir, "chromophores", "output_orca", "single")
    os.makedirs(orca_inp_dir)
    os.makedirs(orca_out_dir)

    try:
        shutil.rmtree(output_dir)
    except OSError:
        pass
    asset_name = os.path.join(TEST_ROOT, "assets", "orca_retrials", "000000.inp"),
    shutil.copy(asset_name, orca_inp_dir)
    file_name = os.path.join(output_dir, os.path.split(asset_name)[1])


    # Create dummy inputs for rerun_fails
    # First the failed_chromo_files dictionary that keeps track of the fail count
    failed_chromo_files = {file_name: fail_count}
    # Then, the parameter dictionary (used to get the orca output directory and the
    # proc_IDs
    parameter_dict = {
        "proc_IDs": [0],
        "output_orca_directory": output_dir
        "output_morphology_directory": output_dir
    }
    # Finally a blank chromophore class (used only to output AAIDs after 18 fails)
    chromophore_list = [dummy_chromophore()]

    ti.rerun_fails(failed_chromo_files, parameter_dict, chromophore_list)


# ---==============================================---
# ---================= Run Tests ==================---
# ---==============================================---


class TestCompareOutputs(TestCommand):
    def test_check_output_exists(self, run_simulation):
        self.confirm_file_exists(
            os.path.join(TEST_ROOT, "output_OR", "chromophores", "output_orca", "single", "000000.inp")
        )


def teardown_module():
    shutil.rmtree(os.path.join(TEST_ROOT, "output_OR"))


if __name__ == "__main__":

    class parameters:
        def __init__(self, param):
            self.param = param

    run_simulation(parameters(""))
