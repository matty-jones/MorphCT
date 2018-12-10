import csv
import os
import pytest
import shutil
import subprocess
from morphct import run_MorphCT
from morphct.definitions import TEST_ROOT
from testing_tools import TestCommand
from morphct.code import helper_functions as hf


@pytest.fixture(scope="module", params=["in", "out"])
def run_simulation(request):
    # ---==============================================---
    # ---=============== Setup Prereqs ================---
    # ---==============================================---
    orca_type = request.param
    orca_file_type = orca_type
    if orca_file_type == "in":
        orca_file_type = "inp"

    output_dir = os.path.join(TEST_ROOT, "output_O2X")

    try:
        shutil.rmtree(output_dir)
    except OSError:
        pass
    os.makedirs(output_dir)
    shutil.copy(
        os.path.join(
            TEST_ROOT,
            "assets",
            "donor_polymer",
            "EZ",
            "".join([orca_type, "put_orca"]),
            "single",
            "".join(["00002.", orca_file_type]),
        ),
        os.path.join(output_dir, "".join(["00002.", orca_file_type])),
    )

    command = [
        "orca2xyz",
        os.path.join(output_dir, "".join(["00002.", orca_file_type])),
    ]
    subprocess.Popen(command).communicate()


# ---==============================================---
# ---================= Run Tests ==================---
# ---==============================================---


class TestCompareOutputs(TestCommand):
    def test_check_output_exists(self, run_simulation):
        self.confirm_file_exists(os.path.join(TEST_ROOT, "output_O2X", "00002.xyz"))


def teardown_module():
    shutil.rmtree(os.path.join(TEST_ROOT, "output_O2X"))


if __name__ == "__main__":

    class parameters:
        def __init__(self, param):
            self.param = param

    run_simulation(parameters("in"))
