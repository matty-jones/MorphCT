from morphct.definitions import TEST_ROOT
from morphct.code import helper_functions as hf
from comparisons import TestCommand
from morphct import run_MorphCT
import os
import shutil
import sys


# ---==============================================---
# ---======== Directory and File Structure ========---
# ---==============================================---

input_morph_dir =   TEST_ROOT + '/assets'
output_morph_dir =  TEST_ROOT + '/temp'
input_device_dir =  TEST_ROOT + '/assets'
output_device_dir = TEST_ROOT + '/temp'

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

execute_fine_graining = True                 # Requires: None
execute_molecular_dynamics = False            # Requires: fine_graining
execute_obtain_chromophores = False           # Requires: Atomistic morphology, or molecular_dynamics
execute_zindo = False                         # Requires: obtain_chromophores
execute_calculate_transfer_integrals = False  # Requires: execute_zindo
execute_calculate_mobility = False            # Requires: calculate_transfer_integrals
execute_device_simulation = False              # Requires: calculate_transfer_integrals for all device_components

# ---==============================================---
# ---========== Fine Graining Parameters ==========---
# ---==============================================---

CG_to_template_dirs = {
    'A': TEST_ROOT + '/assets',
    'B': TEST_ROOT + '/assets',
    'C': TEST_ROOT + '/assets',
}
CG_to_template_files = {
    'A': 'P3HT_template.xml',
    'B': 'P3HT_template.xml',
    'C': 'P3HT_template.xml',
}
CG_to_template_force_fields = {
    'A': 'test_FF.xml',
    'B': 'test_FF.xml',
    'C': 'test_FF.xml',
}
CG_to_template_AAIDs = {
    'A': [0, 1, 2, 3, 4, 24],
    'B': [5, 6, 7, 18, 19, 20, 21, 22, 23],
    'C': [8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
}
CG_to_template_bonds = {
    'bondB': ['C2-C3', 2, 5],
    'bondC': ['C5-C6', 7, 8],
}
rigid_body_sites = {
    'A': [0, 1, 2, 3, 4],
}
additional_constraints = [
    ['C1-C10', 3, 25],
    ['C1-C10-C9', 3, 25, 26],
    ['C1-C10-S1', 3, 25, 29],
    ['S1-C1-C10', 4, 3, 25],
    ['C2-C1-C10', 2, 3, 25],
    ['C1-C10-C9-C2', 3, 25, 26, 27],
    ['C1-C10-S1-C1', 3, 25, 29, 28],
    ['S1-C1-C10-S1', 4, 3, 25, 29],
    ['C2-C1-C10-S1', 2, 3, 25, 29],
]
molecule_terminating_connections = {
    'C1': [[2, 1]],
    'C10': [[2, 1]]
}

# ---==============================================---
# ---================= Begin run ==================---
# ---==============================================---

parameter_file = __file__
proc_IDs = hf.get_CPU_cores()
parameter_names = [i for i in dir() if (not i.startswith('__')) and (not i.startswith('@'))\
                  and (i not in ['run_MorphCT', 'helper_functions', 'hf', 'os', 'shutil', 'TestCommand'])]
parameters = {}
for name in parameter_names:
    parameters[name] = locals()[name]

sys.stdout = None
run_MorphCT.simulation(**parameters)  # Execute MorphCT using these simulation parameters
sys.stdout = sys.__stdout__

# ---==============================================---
# ---================= Run Tests ==================---
# ---==============================================---

def teardown_module(module):
    sys.stdout = sys.__stdout__
    shutil.rmtree(os.path.join(TEST_ROOT, 'temp'))


class TestCompareOutputs(TestCommand):
    # Load the output pickle
    output_pickle_data = hf.load_pickle(os.path.join(output_morph_dir, os.path.splitext(morphology)[0], 'code',
                                                     morphology.replace('.xml', '.pickle')))
    output_AA_morphology_dict = output_pickle_data[0]
    output_CG_morphology_dict = output_pickle_data[1]
    output_CG_to_AAID_master = output_pickle_data[2]
    output_parameter_dict = output_pickle_data[3]
    output_chromophore_list = output_pickle_data[4]

    # Load the expected pickle
    expected_pickle_data = hf.load_pickle(os.path.join(output_morph_dir, os.path.splitext(morphology)[0], 'code',
                                                       morphology.replace('.xml', '.pickle')))
    expected_AA_morphology_dict = expected_pickle_data[0]
    expected_CG_morphology_dict = expected_pickle_data[1]
    expected_CG_to_AAID_master = expected_pickle_data[2]
    expected_parameter_dict = expected_pickle_data[3]
    expected_chromophore_list = expected_pickle_data[4]

    def test_check_AA_morphology_dict(self):
        self.compare_equal(self.output_AA_morphology_dict, self.expected_AA_morphology_dict)

    def test_check_CG_morphology_dict(self):
        self.compare_equal(self.output_CG_morphology_dict,  self.expected_CG_morphology_dict)

    def test_check_CG_to_AAID_master(self):
        self.compare_equal(self.output_CG_to_AAID_master,  self.expected_CG_to_AAID_master)

    def test_check_parameter_dict(self):
        self.compare_equal(self.output_parameter_dict,  self.expected_parameter_dict)

    def test_check_chromophore_list(self):
        self.compare_equal(self.output_chromophore_list,  self.expected_chromophore_list)
