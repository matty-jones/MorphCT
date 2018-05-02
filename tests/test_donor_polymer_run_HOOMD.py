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
    output_morph_dir = TEST_ROOT + '/output_RH'
    input_device_dir = TEST_ROOT + '/assets/donor_polymer'
    output_device_dir = TEST_ROOT + '/output_RH'

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
    execute_molecular_dynamics = True            # Requires: fine_graining
    execute_obtain_chromophores = False           # Requires: Atomistic morphology, or molecular_dynamics
    execute_zindo = False                         # Requires: obtain_chromophores
    execute_calculate_transfer_integrals = False  # Requires: execute_zindo
    execute_calculate_mobility = False            # Requires: calculate_transfer_integrals
    execute_device_simulation = False              # Requires: calculate_transfer_integrals for all device_components

    # ---==============================================---
    # ---=========== Forcefield Parameters ============---
    # ---==============================================---
    CG_to_template_force_fields = {
        'A': 'test_FF.xml',
        'B': 'test_FF.xml',
        'C': 'test_FF.xml',
    }
    pair_R_cut = 10.0
    pair_dpd_gamma_val = 0.0

    # ---==============================================---
    # ---===== Molecular Dynamics Phase Parameters ====---
    # ---==============================================---
    number_of_phases = 4
    temperatures = [1.0]
    taus = [1.0]
    pair_types = ['none', 'dpd', 'lj', 'lj', 'lj', 'lj', 'lj', 'lj']
    bond_types = ['harmonic']
    angle_types = ['harmonic']
    dihedral_types = ['opls']
    integration_targets = ['all']
    timesteps = [1E-3, 1E-3, 1E-7, 1E-4]
    durations = [1E5, 1E5, 1E4, 1E5]
    termination_conditions = ['ke_min', 'max_t', 'max_t', 'max_t']
    group_anchorings = ['all', 'all', 'all', 'none']
    dcd_file_write = True
    dcd_file_dumpsteps = [0]

    # ---==============================================---
    # ---================= Begin run ==================---
    # ---==============================================---

    parameter_file = os.path.realpath(__file__)
    proc_IDs = hf.get_CPU_cores()
    parameter_names = [i for i in dir() if (not i.startswith('__')) and (not i.startswith('@'))\
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
    os.makedirs(os.path.join(output_morph_dir, 'donor_polymer', 'code'))
    shutil.copy(os.path.join(TEST_ROOT, 'assets', 'donor_polymer', 'FG',
                             'donor_polymer_post_fine_graining.pickle'),
                os.path.join(output_morph_dir, 'donor_polymer', 'code', 'donor_polymer.pickle'))

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
    expected_pickle_data = hf.load_pickle(os.path.join(input_morph_dir, 'RH',
                                                       'donor_polymer_post_run_HOOMD.pickle'))
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

    def test_check_morphology_phase0_created(self, run_simulation):
        output_morph_dir = run_simulation['output_parameter_dict']['output_morph_dir']
        morphology = run_simulation['output_parameter_dict']['morphology']
        morph_dir = os.path.join(output_morph_dir, os.path.splitext(morphology)[0], 'morphology')
        self.confirm_file_exists(os.path.join(morph_dir, 'phase0_' + morphology))

    def test_check_morphology_phase1_created(self, run_simulation):
        output_morph_dir = run_simulation['output_parameter_dict']['output_morph_dir']
        morphology = run_simulation['output_parameter_dict']['morphology']
        morph_dir = os.path.join(output_morph_dir, os.path.splitext(morphology)[0], 'morphology')
        self.confirm_file_exists(os.path.join(morph_dir, 'phase1_' + morphology))
        self.confirm_file_exists(os.path.join(morph_dir, 'phase1_' + morphology.replace('xml', 'dcd')))

    def test_check_morphology_phase2_created(self, run_simulation):
        output_morph_dir = run_simulation['output_parameter_dict']['output_morph_dir']
        morphology = run_simulation['output_parameter_dict']['morphology']
        morph_dir = os.path.join(output_morph_dir, os.path.splitext(morphology)[0], 'morphology')
        self.confirm_file_exists(os.path.join(morph_dir, 'phase2_' + morphology))
        self.confirm_file_exists(os.path.join(morph_dir, 'phase2_' + morphology.replace('xml', 'dcd')))

    def test_check_morphology_phase3_created(self, run_simulation):
        output_morph_dir = run_simulation['output_parameter_dict']['output_morph_dir']
        morphology = run_simulation['output_parameter_dict']['morphology']
        morph_dir = os.path.join(output_morph_dir, os.path.splitext(morphology)[0], 'morphology')
        self.confirm_file_exists(os.path.join(morph_dir, 'phase3_' + morphology))
        self.confirm_file_exists(os.path.join(morph_dir, 'phase3_' + morphology.replace('xml', 'dcd')))

    def test_check_morphology_phase4_created(self, run_simulation):
        output_morph_dir = run_simulation['output_parameter_dict']['output_morph_dir']
        morphology = run_simulation['output_parameter_dict']['morphology']
        morph_dir = os.path.join(output_morph_dir, os.path.splitext(morphology)[0], 'morphology')
        self.confirm_file_exists(os.path.join(morph_dir, 'phase4_' + morphology))
        self.confirm_file_exists(os.path.join(morph_dir, 'phase4_' + morphology.replace('xml', 'dcd')))

    def test_check_morphology_final_created(self, run_simulation):
        output_morph_dir = run_simulation['output_parameter_dict']['output_morph_dir']
        morphology = run_simulation['output_parameter_dict']['morphology']
        morph_dir = os.path.join(output_morph_dir, os.path.splitext(morphology)[0], 'morphology')
        self.confirm_file_exists(os.path.join(morph_dir, 'final_' + morphology))
        self.confirm_file_exists(os.path.join(morph_dir, 'energies_' + morphology.replace('xml', 'log')))


def teardown_module():
    shutil.rmtree(os.path.join(TEST_ROOT, 'output_RH'))


if __name__ == "__main__":
    run_simulation()
