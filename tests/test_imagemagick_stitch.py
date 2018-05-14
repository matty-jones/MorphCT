import csv
import os
import pytest
import shutil
import subprocess
from morphct import run_MorphCT
from morphct.definitions import TEST_ROOT
from testing_tools import TestCommand
from morphct.code import helper_functions as hf


@pytest.fixture(scope='module',
                params=['', '-t test', '-d 1x', '-d 2x2', '-d x1']
                )
def run_simulation(request):
    flags = request.param

    # ---==============================================---
    # ---=============== Setup Prereqs ================---
    # ---==============================================---

    output_dir = os.path.join(TEST_ROOT, 'output_IM')

    try:
        shutil.rmtree(output_dir)
    except OSError:
        pass
    #os.makedirs(os.path.join(output_dir, 'figures'))
    shutil.copytree(os.path.join(TEST_ROOT, 'assets', 'imagemagick_stitch'),
                    os.path.join(output_dir, 'figures'))
    command = [component for component in ['createMontage'] + flags.split()
               + [output_dir]]
    print("Executing command", command)
    subprocess.Popen(command).communicate()
    return flags


# ---==============================================---
# ---================= Run Tests ==================---
# ---==============================================---


class TestCompareOutputs(TestCommand):
    def test_check_output_exists(self, run_simulation):
        if run_simulation[:2] == '-t':
            file_name = ''.join([run_simulation[3:], '.png'])
        else:
            file_name = 'output_IM.png'
        self.confirm_file_exists(os.path.join(
            TEST_ROOT, 'output_IM', file_name))


def teardown_module():
    shutil.rmtree(os.path.join(TEST_ROOT, 'output_IM'))


if __name__ == "__main__":
    class parameters:
        def __init__(self, param):
            self.param = param

    run_simulation(parameters('-t test'))
