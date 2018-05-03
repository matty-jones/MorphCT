from morphct.definitions import TEST_ROOT
from morphct.code import helper_functions as hf
from morphct import run_MorphCT
from testing_tools import TestCommand
import os
import shutil
import sys
import pytest

@pytest.fixture(scope='module')
def run_simulation():
    # ---==============================================---
    # ---======== Directory and File Structure ========---
    # ---==============================================---

    input_morph_dir = TEST_ROOT + '/assets/donor_polymer'
    output_morph_dir = TEST_ROOT + '/output_TI'
    input_device_dir = TEST_ROOT + '/assets/donor_polymer'
    output_device_dir = TEST_ROOT + '/output_TI'

    # ---==============================================---
    # ---========== Input Morphology Details ==========---
    # ---==============================================---

    morphology = 'donor_polymer.xml'
    input_sigma = 3.0
    device_morphology = None
    device_components = {
    }
    overwrite_current_data = True

    # ---==============================================---
    # ---============= Execution Modules ==============---
    # ---==============================================---

    execute_fine_graining = False                 # Requires: None
    execute_molecular_dynamics = False            # Requires: fine_graining
    execute_obtain_chromophores = False           # Requires: Atomistic morphology, or molecular_dynamics
    execute_zindo = False                         # Requires: obtain_chromophores
    execute_calculate_transfer_integrals = True  # Requires: execute_zindo
    execute_calculate_mobility = False            # Requires: calculate_transfer_integrals
    execute_device_simulation = False              # Requires: calculate_transfer_integrals for all device_components


    # ---==============================================---
    # ---================= Begin run ==================---
    # ---==============================================---

    parameter_file = os.path.realpath(__file__)
    proc_IDs = hf.get_CPU_cores()
    parameter_names = [i for i in dir() if (not i.startswith('__')) and (not i.startswith('@'))\
                       and (not i.startswith('Test')) and (not i.startswith('test'))\
                       and (i not in ['run_MorphCT', 'helper_functions', 'hf', 'os', 'shutil', 'TestCommand',
                                      'TEST_ROOT', 'setup_module', 'teardown_module', 'testing_tools', 'sys',
                                      'pytest'])]
    parameters = {}
    for name in parameter_names:
        parameters[name] = locals()[name]

    # ---==============================================---
    # ---=============== Setup Prereqs ================---
    # ---==============================================---
    try:
        shutil.rmtree(output_morph_dir)
    except OSError:
        pass
    os.makedirs(os.path.join(output_morph_dir, os.path.splitext(morphology)[0], 'code'))
    shutil.copy(os.path.join(TEST_ROOT, 'assets', os.path.splitext(morphology)[0], 'EZ',
                             morphology.replace('.xml', '_post_execute_ZINDO.pickle')),
                os.path.join(output_morph_dir, os.path.splitext(morphology)[0], 'code',
                             morphology.replace('.xml', '.pickle')))
    shutil.copytree(os.path.join(TEST_ROOT, 'assets', os.path.splitext(morphology)[0], 'EZ',
                                 'input_orca'),
                os.path.join(output_morph_dir, os.path.splitext(morphology)[0], 'chromophores',
                             'input_orca'))
    shutil.copytree(os.path.join(TEST_ROOT, 'assets', os.path.splitext(morphology)[0], 'EZ',
                                 'output_orca'),
                os.path.join(output_morph_dir, os.path.splitext(morphology)[0], 'chromophores',
                             'output_orca'))

    run_MorphCT.simulation(**parameters)  # Execute MorphCT using these simulation parameters
    # The output dictionary from this fixing
    fix_dict = {}

    # Load the output pickle
    output_pickle_data = hf.load_pickle(os.path.join(output_morph_dir, os.path.splitext(morphology)[0],
                                                          'code', morphology.replace('.xml', '.pickle')))
    fix_dict['output_AA_morphology_dict'] = output_pickle_data[0]
    fix_dict['output_CG_morphology_dict'] = output_pickle_data[1]
    fix_dict['output_CG_to_AAID_master'] = output_pickle_data[2]
    fix_dict['output_parameter_dict'] = output_pickle_data[3]
    fix_dict['output_chromophore_list'] = output_pickle_data[4]

    # Load the expected pickle
    expected_pickle_data = hf.load_pickle(
        os.path.join(input_morph_dir, 'TI',
                     morphology.replace('.xml', '_post_calculate_transfer_integrals.pickle')))
    fix_dict['expected_AA_morphology_dict'] = expected_pickle_data[0]
    fix_dict['expected_CG_morphology_dict'] = expected_pickle_data[1]
    fix_dict['expected_CG_to_AAID_master'] = expected_pickle_data[2]
    fix_dict['expected_parameter_dict'] = expected_pickle_data[3]
    fix_dict['expected_chromophore_list'] = expected_pickle_data[4]
    return fix_dict

# ---==============================================---
# ---================= Run Tests ==================---
# ---==============================================---


class TestCompareOutputs(TestCommand):
    def test_check_AA_morphology_dict(self, run_simulation):
        self.compare_equal(run_simulation['output_AA_morphology_dict'],
                           run_simulation['expected_AA_morphology_dict'])

    def test_check_CG_morphology_dict(self, run_simulation):
        self.compare_equal(run_simulation['output_CG_morphology_dict'],
                           run_simulation['expected_CG_morphology_dict'])

    def test_check_CG_to_AAID_master(self, run_simulation):
        self.compare_equal(run_simulation['output_CG_to_AAID_master'],
                           run_simulation['expected_CG_to_AAID_master'])

    def test_check_parameter_dict(self, run_simulation):
        # Pop the system-dependent keys, such as the input and output dirs since this will
        # always be system-dependent
        output_pars = {}
        expected_pars = {}
        for key in run_simulation['expected_parameter_dict']:
            if key  in ['parameter_file', 'output_morph_dir', 'CG_to_template_dirs', 'output_morphology_directory',
                        'input_device_dir', 'input_morphology_file', 'output_device_dir', 'input_morph_dir']:
                continue
            output_pars = run_simulation['output_parameter_dict'][key]
            expected_pars = run_simulation['expected_parameter_dict'][key]
        self.compare_equal(output_pars, expected_pars)

    def test_check_chromophore_list(self, run_simulation):
        self.compare_equal(run_simulation['output_chromophore_list'],
                           run_simulation['expected_chromophore_list'])


# TODO: Tests for failed singles and failed pairs


def teardown_module():
    shutil.rmtree(os.path.join(TEST_ROOT, 'output_TI'))


if __name__ == "__main__":
    run_simulation()
