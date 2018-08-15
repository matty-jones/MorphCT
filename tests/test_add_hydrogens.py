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
        " ".join(["-m", "DBP", os.path.join(TEST_ROOT, "output_AH", "DBP.xml")]),
        " ".join(["-m", "P3HT", os.path.join(TEST_ROOT, "output_AH", "P3HT.xml")]),
        " ".join(
            ["-m", "P3HT/PCBM", os.path.join(TEST_ROOT, "output_AH", "PC60BM.xml")]
        ),
        " ".join(
            ["-m", "P3HT/PCBM", os.path.join(TEST_ROOT, "output_AH", "PC70BM.xml")]
        ),
        " ".join(["-m", "PERYLENE", os.path.join(TEST_ROOT, "output_AH", "PE.xml")]),
        " ".join(["-m", "PERYLENE", os.path.join(TEST_ROOT, "output_AH", "PT.xml")]),
        " ".join(["-m", "BDT-TPD", os.path.join(TEST_ROOT, "output_AH", "BDTTPD.xml")]),
    ],
)
def run_simulation(request):
    flags = request.param

    # ---==============================================---
    # ---=============== Setup Prereqs ================---
    # ---==============================================---

    output_dir = os.path.join(TEST_ROOT, "output_AH")

    try:
        shutil.rmtree(output_dir)
    except OSError:
        pass
    os.makedirs(os.path.join(output_dir))
    molecule_file = os.path.split(flags)[1]
    shutil.copy(
        os.path.join(TEST_ROOT, "assets", "add_hydrogens", molecule_file),
        os.path.join(output_dir),
    )
    command = [component for component in ["addHydrogens"] + flags.split()]
    print("Executing command", command)
    subprocess.Popen(command).communicate()
    return flags


# ---==============================================---
# ---================= Run Tests ==================---
# ---==============================================---


class TestCompareOutputs(TestCommand):
    def test_check_output_exists(self, run_simulation):
        output_file = os.path.split(run_simulation)[1].replace(".xml", "_AA.xml")
        response_location = os.path.join(TEST_ROOT, "output_AH", output_file)
        self.confirm_file_exists(response_location)

    def test_check_output_identical(self, run_simulation):
        output_file = os.path.split(run_simulation)[1].replace(".xml", "_AA.xml")
        response_location = os.path.join(TEST_ROOT, "output_AH", output_file)
        expected_location = os.path.join(
            TEST_ROOT, "assets", "add_hydrogens", output_file
        )
        response_morph = hf.load_morphology_xml(response_location)
        expected_morph = hf.load_morphology_xml(expected_location)
        self.compare_equal(expected_morph, response=response_morph)


def teardown_module():
    shutil.rmtree(os.path.join(TEST_ROOT, "output_AH"))


if __name__ == "__main__":

    class parameters:
        def __init__(self, param):
            self.param = param

    run_simulation(parameters("-m BDT-TPD output_AH/BDTTPD.xml"))
    # run_simulation(parameters('-m DBP output_AH/DBP.xml'))
    # run_simulation(parameters('-m P3HT output_AH/P3HT.xml'))
    # run_simulation(parameters('-m P3HT/PCBM output_AH/PC60BM.xml'))
    # run_simulation(parameters('-m P3HT/PCBM output_AH/PC70BM.xml'))
    # run_simulation(parameters('-m PERYLENE output_AH/PE.xml'))
    # run_simulation(parameters('-m PERYLENE output_AH/PT.xml'))
