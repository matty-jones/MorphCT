import copy
import os
import pytest
import shutil
from morphct import run_MorphCT
from morphct.definitions import TEST_ROOT
from testing_tools import TestCommand
from morphct.code import helper_functions as hf


@pytest.fixture(scope="module")
def run_simulation():
    # ---==============================================---
    # ---======== Directory and File Structure ========---
    # ---==============================================---

    input_morph_dir = os.path.join(TEST_ROOT, "assets", "donor_polymer")
    output_morph_dir = os.path.join(TEST_ROOT, "output_FG")
    output_orca_dir = None
    input_device_dir = os.path.join(TEST_ROOT, "assets", "donor_polymer")
    output_device_dir = os.path.join(TEST_ROOT, "output_FG")

    # ---==============================================---
    # ---========== Input Morphology Details ==========---
    # ---==============================================---

    morphology = "donor_polymer.xml"
    input_sigma = 3.0
    device_morphology = None
    device_components = {}
    overwrite_current_data = True
    random_seed_override = 929292929

    # ---==============================================---
    # ---============= Execution Modules ==============---
    # ---==============================================---

    execute_fine_graining = True  # Requires: None
    execute_molecular_dynamics = False  # Requires: fine_graining
    execute_obtain_chromophores = (
        False
    )  # Requires: Atomistic morphology, or molecular_dynamics
    execute_ZINDO = False  # Requires: obtain_chromophores
    execute_calculate_transfer_integrals = False  # Requires: execute_ZINDO
    execute_calculate_mobility = False  # Requires: calculate_transfer_integrals
    execute_device_simulation = (
        False
    )  # Requires: calculate_transfer_integrals for all device_components

    # ---==============================================---
    # ---========== Fine Graining Parameters ==========---
    # ---==============================================---

    CG_to_template_dirs = {
        "A": os.path.join(TEST_ROOT, "assets", "donor_polymer"),
        "B": os.path.join(TEST_ROOT, "assets", "donor_polymer"),
        "C": os.path.join(TEST_ROOT, "assets", "donor_polymer"),
    }
    CG_to_template_files = {
        "A": "P3HT_template.xml",
        "B": "P3HT_template.xml",
        "C": "P3HT_template.xml",
    }
    CG_to_template_force_fields = {
        "A": "test_FF.xml",
        "B": "test_FF.xml",
        "C": "test_FF.xml",
    }
    CG_to_template_AAIDs = {
        "A": [0, 1, 2, 3, 4, 24],
        "B": [5, 6, 7, 18, 19, 20, 21, 22, 23],
        "C": [8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
    }
    CG_to_template_bonds = {"bondB": ["C2-C3", 2, 5], "bondC": ["C5-C6", 7, 8]}
    rigid_body_sites = {"A": [0, 1, 2, 3, 4]}
    additional_constraints = [
        ["C1-C10", 3, 25],
        ["C1-C10-C9", 3, 25, 26],
        ["C1-C10-S1", 3, 25, 29],
        ["S1-C1-C10", 4, 3, 25],
        ["C2-C1-C10", 2, 3, 25],
        ["C1-C10-C9-C2", 3, 25, 26, 27],
        ["C1-C10-S1-C1", 3, 25, 29, 28],
        ["S1-C1-C10-S1", 4, 3, 25, 29],
        ["C2-C1-C10-S1", 2, 3, 25, 29],
    ]
    molecule_terminating_connections = {"C1": [[2, 1]], "C10": [[2, 1]]}

    # ---==============================================---
    # ---================= Begin run ==================---
    # ---==============================================---

    parameter_file = os.path.realpath(__file__)
    proc_IDs = hf.get_CPU_cores()
    parameter_names = [
        i
        for i in dir()
        if (not i.startswith("__"))
        and (not i.startswith("@"))
        and (not i.startswith("Test"))
        and (not i.startswith("test"))
        and (
            i
            not in [
                "run_MorphCT",
                "helper_functions",
                "hf",
                "os",
                "shutil",
                "TestCommand",
                "TEST_ROOT",
                "setup_module",
                "teardown_module",
                "testing_tools",
                "sys",
                "pytest",
            ]
        )
    ]
    parameters = {}
    for name in parameter_names:
        parameters[name] = locals()[name]
    run_MorphCT.simulation(
        **parameters
    )  # Execute MorphCT using these simulation parameters
    # The output dictionary from this fixing
    fix_dict = {}

    # Load the output pickle
    output_pickle_data = hf.load_pickle(
        os.path.join(
            output_morph_dir,
            os.path.splitext(morphology)[0],
            "code",
            morphology.replace(".xml", ".pickle"),
        )
    )
    fix_dict["output_AA_morphology_dict"] = output_pickle_data[0]
    fix_dict["output_CG_morphology_dict"] = output_pickle_data[1]
    fix_dict["output_CG_to_AAID_master"] = output_pickle_data[2]
    fix_dict["output_parameter_dict"] = output_pickle_data[3]
    fix_dict["output_chromophore_list"] = output_pickle_data[4]

    # Load the expected pickle
    expected_pickle_data = hf.load_pickle(
        os.path.join(
            input_morph_dir,
            "FG",
            morphology.replace(".xml", "_post_fine_graining.pickle"),
        )
    )
    fix_dict["expected_AA_morphology_dict"] = expected_pickle_data[0]
    fix_dict["expected_CG_morphology_dict"] = expected_pickle_data[1]
    fix_dict["expected_CG_to_AAID_master"] = expected_pickle_data[2]
    fix_dict["expected_parameter_dict"] = expected_pickle_data[3]
    fix_dict["expected_chromophore_list"] = expected_pickle_data[4]
    return fix_dict


# ---==============================================---
# ---================= Run Tests ==================---
# ---==============================================---


class TestCompareOutputs(TestCommand):
    def test_check_AA_morphology_dict(self, run_simulation):
        self.compare_equal(
            run_simulation["expected_AA_morphology_dict"],
            response=run_simulation["output_AA_morphology_dict"],
        )

    def testck_CG_morphology_dict(self, run_simulation):
        self.compare_equal(
            run_simulation["expected_CG_morphology_dict"],
            response=run_simulation["output_CG_morphology_dict"],
        )

    def test_check_CG_to_AAID_master(self, run_simulation):
        self.compare_equal(
            run_simulation["expected_CG_to_AAID_master"],
            response=run_simulation["output_CG_to_AAID_master"],
        )

    def test_check_parameter_dict(self, run_simulation):
        # Pop the system-dependent keys, such as the input and output dirs since this will
        # always be system-dependent
        expected_pars = copy.deepcopy(run_simulation["expected_parameter_dict"])
        output_pars = copy.deepcopy(run_simulation["output_parameter_dict"])
        for key in [
            "parameter_file",
            "output_morph_dir",
            "CG_to_template_dirs",
            "output_morphology_directory",
            "input_device_dir",
            "input_morphology_file",
            "output_device_dir",
            "input_morph_dir",
            "input_orca_dir",
            "output_orca_dir",
            "input_device_file",
            "output_device_directory",
            "output_orca_directory",
        ]:
            try:
                expected_pars.pop(key)
            except KeyError:
                pass
            try:
                output_pars.pop(key)
            except KeyError:
                pass
        self.compare_equal(expected_pars, response=output_pars)

    def test_check_chromophore_list(self, run_simulation):
        self.compare_equal(
            run_simulation["expected_chromophore_list"],
            response=run_simulation["output_chromophore_list"],
        )

    def test_check_code_copied(self, run_simulation):
        output_morph_dir = run_simulation["output_parameter_dict"]["output_morph_dir"]
        morphology = run_simulation["output_parameter_dict"]["morphology"]
        code_dir = os.path.join(
            output_morph_dir, os.path.splitext(morphology)[0], "code"
        )
        self.confirm_file_exists(os.path.join(code_dir, "fine_grainer.py"))
        self.confirm_file_exists(os.path.join(code_dir, "run_HOOMD.py"))
        self.confirm_file_exists(os.path.join(code_dir, "obtain_chromophores.py"))
        self.confirm_file_exists(os.path.join(code_dir, "execute_ZINDO.py"))
        self.confirm_file_exists(os.path.join(code_dir, "transfer_integrals.py"))
        self.confirm_file_exists(os.path.join(code_dir, "mobility_KMC.py"))
        self.confirm_file_exists(os.path.join(code_dir, "device_KMC.py"))
        self.confirm_file_exists(
            os.path.join(code_dir, "single_core_run_device_KMC.py")
        )
        self.confirm_file_exists(os.path.join(code_dir, "single_core_run_mob_KMC.py"))
        self.confirm_file_exists(os.path.join(code_dir, "single_core_run_orca.py"))

    def test_check_par_copied(self, run_simulation):
        output_morph_dir = run_simulation["output_parameter_dict"]["output_morph_dir"]
        morphology = run_simulation["output_parameter_dict"]["morphology"]
        code_dir = os.path.join(
            output_morph_dir, os.path.splitext(morphology)[0], "code"
        )
        self.confirm_file_exists(os.path.join(code_dir, __file__))

    def test_check_input_morph_copied(self, run_simulation):
        output_morph_dir = run_simulation["output_parameter_dict"]["output_morph_dir"]
        morphology = run_simulation["output_parameter_dict"]["morphology"]
        code_dir = os.path.join(
            output_morph_dir, os.path.splitext(morphology)[0], "code"
        )
        self.confirm_file_exists(os.path.join(code_dir, "input.xml"))

    def test_check_morphology_created(self, run_simulation):
        output_morph_dir = run_simulation["output_parameter_dict"]["output_morph_dir"]
        morphology = run_simulation["output_parameter_dict"]["morphology"]
        morph_dir = os.path.join(
            output_morph_dir, os.path.splitext(morphology)[0], "morphology"
        )
        self.confirm_file_exists(os.path.join(morph_dir, morphology))

    def test_check_morphology(self, run_simulation):
        output_morph_dir = run_simulation["output_parameter_dict"]["output_morph_dir"]
        input_morph_dir = run_simulation["output_parameter_dict"]["input_morph_dir"]
        morphology = run_simulation["output_parameter_dict"]["morphology"]
        morph_dir = os.path.join(
            output_morph_dir, os.path.splitext(morphology)[0], "morphology"
        )
        output_morphology = hf.load_morphology_xml(os.path.join(morph_dir, morphology))
        expected_morphology = hf.load_morphology_xml(
            os.path.join(
                input_morph_dir,
                "FG",
                morphology.replace(".xml", "_post_fine_graining.xml"),
            )
        )
        self.compare_equal(expected_morphology, response=output_morphology)


def teardown_module():
    shutil.rmtree(os.path.join(TEST_ROOT, "output_FG"))


if __name__ == "__main__":
    run_simulation()
