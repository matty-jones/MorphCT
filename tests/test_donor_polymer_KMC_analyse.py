import os
import pytest
import shutil
import subprocess
from morphct import run_MorphCT
from morphct.definitions import TEST_ROOT
from testing_tools import TestCommand
from morphct.code import helper_functions as hf


@pytest.fixture(scope='module',
                params=['', '-tp']
                )
def run_simulation(request):
    flags = request.param

    # ---==============================================---
    # ---=============== Setup Prereqs ================---
    # ---==============================================---

    output_dir = 'output_KMCA'

    try:
        shutil.rmtree(output_dir)
    except OSError:
        pass
    os.makedirs(os.path.join(output_dir, 'code'))
    os.makedirs(os.path.join(output_dir, 'figures'))
    shutil.copy(os.path.join(TEST_ROOT, 'assets', 'donor_polymer', 'MKMC',
                             'donor_polymer_post_calculate_mobility.pickle'),
                os.path.join(output_dir, 'code', ''.join([output_dir,
                                                          '.pickle'])))
    shutil.copytree(os.path.join(TEST_ROOT, 'assets', 'donor_polymer', 'MKMC',
                                 'KMC'),
                    os.path.join(output_dir, 'KMC'))
    command = [component for component in ['KMCAnalyse', flags, output_dir] if
               len(component) > 0]
    print("Executing command", command)
    subprocess.Popen(command).communicate()
    exit()
    return 0


# ---==============================================---
# ---================= Run Tests ==================---
# ---==============================================---


class TestCompareOutputs(TestCommand):
    def test_check_AA_morphology_dict(self, run_simulation):
        return True


def teardown_module():
    shutil.rmtree(os.path.join(TEST_ROOT, 'output_KMCA'))


if __name__ == "__main__":
    class parameters:
        def __init__(self, param):
            self.param = param

    run_simulation(parameters('-tp'))
