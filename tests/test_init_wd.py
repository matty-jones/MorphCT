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

    command = ["morphctInit"]
    subprocess.Popen(command).communicate()


# ---==============================================---
# ---================= Run Tests ==================---
# ---==============================================---


class TestCompareOutputs(TestCommand):
    def test_check_dirs_exist(self, run_simulation):
        os.path.isdir(os.path.join(os.getcwd(), "inputs"))
        os.path.isdir(os.path.join(os.getcwd(), "outputs"))

    def test_check_files_exist(self, run_simulation):
        self.confirm_file_exists(os.path.join(os.getcwd(), "submit.sh"))
        self.confirm_file_exists(os.path.join(os.getcwd(), "submit.py"))
        self.confirm_file_exists(os.path.join(os.getcwd(), "par.py"))


def teardown_module():
    shutil.rmtree(os.path.join(os.getcwd(), "inputs"))
    shutil.rmtree(os.path.join(os.getcwd(), "outputs"))
    os.remove(os.path.join(os.getcwd(), "submit.sh"))
    os.remove(os.path.join(os.getcwd(), "submit.py"))
    os.remove(os.path.join(os.getcwd(), "par.py"))


if __name__ == "__main__":

    class parameters:
        def __init__(self, param):
            self.param = param

    run_simulation(parameters(""))
