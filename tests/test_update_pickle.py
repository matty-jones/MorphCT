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
        "output_UP/MCT1.0_pickle",
        "-p output_UP/MCT2.0_pickle/MCT2.0_pickle.pickle",
        "output_UP/MCT3.0_pickle",
    ],
)
def run_simulation(request):
    # ---==============================================---
    # ---=============== Setup Prereqs ================---
    # ---==============================================---
    flags = request.param

    output_dir = os.path.join(TEST_ROOT, "output_UP")

    try:
        shutil.rmtree(output_dir)
    except OSError:
        pass
    # os.makedirs(os.path.join(output_dir, 'figures'))
    shutil.copytree(
        os.path.join(TEST_ROOT, "assets", "update_pickle"), os.path.join(output_dir)
    )

    command = ["updatePickle"] + flags.split()
    subprocess.Popen(command).communicate()
    return os.path.split(flags.split()[-1])


# ---==============================================---
# ---================= Run Tests ==================---
# ---==============================================---


class TestCompareOutputs(TestCommand):
    def test_check_output_exists(self, run_simulation):
        (directory, pickle_file) = run_simulation
        self.confirm_file_exists(
            os.path.join(
                TEST_ROOT, directory, pickle_file.replace(".pickle", ".bak_pickle")
            )
        )


def teardown_module():
    shutil.rmtree(os.path.join(TEST_ROOT, "output_UP"))


if __name__ == "__main__":

    class parameters:
        def __init__(self, param):
            self.param = param
    run_simulation(parameters("output_UP/MCT1.0_pickle"))
