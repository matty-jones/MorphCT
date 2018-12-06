import csv
import os
import pytest
import shutil
import subprocess
from morphct import run_MorphCT
from morphct.definitions import TEST_ROOT
from testing_tools import TestCommand
from morphct.code import helper_functions as hf


@pytest.fixture(
    scope="module",
    params=[
        os.path.join(TEST_ROOT, "output_UP/MCT1.0_pickle"),
        " ".join([
            "-p",
            os.path.join(TEST_ROOT, "output_UP/MCT2.0_pickle/MCT2.0_pickle.pickle")
        ]),
        os.path.join(TEST_ROOT, "output_UP/MCT3.0_pickle"),
    ],
)
def run_simulation(request):
    # ---==============================================---
    # ---=============== Setup Prereqs ================---
    # ---==============================================---
    flags = request.param.split()

    output_dir = os.path.join(TEST_ROOT, "output_UP")

    try:
        shutil.rmtree(output_dir)
    except OSError:
        pass
    # os.makedirs(os.path.join(output_dir, 'figures'))
    shutil.copytree(
        os.path.join(TEST_ROOT, "assets", "update_pickle"), os.path.join(output_dir)
    )

    command = ["updatePickle"] + flags
    subprocess.Popen(command).communicate()
    return os.path.split(flags[-1])


# ---==============================================---
# ---================= Run Tests ==================---
# ---==============================================---


class TestCompareOutputs(TestCommand):
    def test_check_outputs_exist(self, run_simulation):
        (directory, pickle_file) = run_simulation
        self.confirm_file_exists(
            os.path.join(
                TEST_ROOT, directory, pickle_file.replace(".pickle", ".bak_pickle")
            )
        )
        if "MCT3.0" in pickle_file:
            # Check for the updated parameter file
            self.confirm_file_exists(
                os.path.join(
                    TEST_ROOT, directory, "MCT3.0_pickle", "code",
                    "MCT3.0_pickle_par.bak_par"
                )
            )

    def test_outputs_import(self, run_simulation):
        (directory, pickle_file) = run_simulation
        if ".pickle" not in pickle_file:
            pickle_file = "".join([pickle_file, "/code/", pickle_file, ".pickle"])
        hf.load_pickle(os.path.join(directory, pickle_file))


def teardown_module():
    shutil.rmtree(os.path.join(TEST_ROOT, "output_UP"))


if __name__ == "__main__":

    class parameters:
        def __init__(self, param):
            self.param = param
    run_simulation(parameters(os.path.join(TEST_ROOT, "output_UP/MCT1.0_pickle")))
