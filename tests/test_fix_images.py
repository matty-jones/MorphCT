import csv
import os
import pytest
import shutil
import subprocess
from morphct import run_MorphCT
from morphct.definitions import TEST_ROOT
from testing_tools import TestCommand
from morphct.code import helper_functions as hf


@pytest.fixture(scope="module")
def run_simulation(request):
    # ---==============================================---
    # ---=============== Setup Prereqs ================---
    # ---==============================================---

    output_dir = os.path.join(TEST_ROOT, "output_FI")

    try:
        shutil.rmtree(output_dir)
    except OSError:
        pass
    # os.makedirs(os.path.join(output_dir, 'figures'))
    shutil.copytree(
        os.path.join(TEST_ROOT, "assets", "fix_images"), os.path.join(output_dir)
    )
    command = ["fixImages", os.path.join(output_dir, "fix_images_input.xml")]
    subprocess.Popen(command).communicate()


# ---==============================================---
# ---================= Run Tests ==================---
# ---==============================================---


class TestCompareOutputs(TestCommand):
    def test_check_output_exists(self, run_simulation):
        self.confirm_file_exists(
            os.path.join(TEST_ROOT, "output_FI", "image_fix_fix_images_input.xml")
        )

    def test_compare_results(self, run_simulation):
        results_dict = hf.load_morphology_xml(
            os.path.join(TEST_ROOT, "output_FI", "image_fix_fix_images_input.xml")
        )
        expected_dict = hf.load_morphology_xml(
            os.path.join(TEST_ROOT, "assets", "fix_images", "fix_images_output.xml")
        )
        self.compare_equal(expected_dict, response=results_dict)


def teardown_module():
    shutil.rmtree(os.path.join(TEST_ROOT, "output_FI"))


if __name__ == "__main__":

    class parameters:
        def __init__(self, param):
            self.param = param

    run_simulation(parameters(""))
