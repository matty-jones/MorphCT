import copy
import os
import pickle
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
    output_morph_dir = os.path.join(TEST_ROOT, "output_MKMC")
    output_orca_dir = None
    input_device_dir = os.path.join(TEST_ROOT, "assets", "donor_polymer")
    output_device_dir = os.path.join(TEST_ROOT, "output_MKMC")

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

    execute_fine_graining = False  # Requires: None
    execute_molecular_dynamics = False  # Requires: fine_graining
    execute_obtain_chromophores = (
        False
    )  # Requires: Atomistic morphology, or molecular_dynamics
    execute_ZINDO = False  # Requires: obtain_chromophores
    execute_calculate_transfer_integrals = False  # Requires: execute_ZINDO
    execute_calculate_mobility = True  # Requires: calculate_transfer_integrals
    execute_device_simulation = (
        False
    )  # Requires: calculate_transfer_integrals for all device_components

    # ---==============================================---
    # ---=== General Kinetic Monte Carlo Parameters ===---
    # ---==============================================---

    # ---=== Universal KMC Parameters ===---
    system_temperature = 290
    use_simple_energetic_penalty = False
    record_carrier_history = True
    use_VRH = True

    # ---=== Mobility Specific KMC Parameters ===---
    number_of_holes_per_simulation_time = 10
    number_of_electrons_per_simulation_time = 0
    hop_limit = 0
    simulation_times = [1.00e-13, 1.00e-12, 1.00e-11]
    combine_KMC_results = True
    use_average_hop_rates = False
    average_intra_hop_rate = 8.07E14
    average_inter_hop_rate = 3.92E14

    # ---==============================================---
    # ---================= Begin run ==================---
    # ---==============================================---

    parameter_file = os.path.realpath(__file__)
    # Force serial running
    proc_IDs = [0]
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
                "pickle",
            ]
        )
    ]
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
    os.makedirs(os.path.join(output_morph_dir, os.path.splitext(morphology)[0], "code"))
    shutil.copy(
        os.path.join(
            TEST_ROOT,
            "assets",
            os.path.splitext(morphology)[0],
            "MKMC",
            morphology.replace(".xml", "_post_calculate_transfer_integrals.pickle"),
        ),
        os.path.join(
            output_morph_dir,
            os.path.splitext(morphology)[0],
            "code",
            morphology.replace(".xml", ".pickle"),
        ),
    )

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
            "MKMC",
            morphology.replace(".xml", "_post_calculate_mobility.pickle"),
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

    def test_check_CG_morphology_dict(self, run_simulation):
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

    def test_check_KMC_outputs_created(self, run_simulation):
        output_morph_dir = run_simulation["output_parameter_dict"]["output_morph_dir"]
        morphology = run_simulation["output_parameter_dict"]["morphology"]
        KMC_dir = os.path.join(output_morph_dir, os.path.splitext(morphology)[0], "KMC")
        self.confirm_file_exists(os.path.join(KMC_dir, "KMC_log_00.log"))
        self.confirm_file_exists(os.path.join(KMC_dir, "KMC_results.pickle"))

    def test_check_KMC_log_identical(self, run_simulation):
        output_morph_dir = run_simulation["output_parameter_dict"]["output_morph_dir"]
        morphology = run_simulation["output_parameter_dict"]["morphology"]
        KMC_dir = os.path.join(output_morph_dir, os.path.splitext(morphology)[0], "KMC")
        with open(os.path.join(KMC_dir, "KMC_log_00.log"), "r") as result_file, open(
            os.path.join(
                TEST_ROOT,
                "assets",
                os.path.splitext(morphology)[0],
                "MKMC",
                "KMC",
                "KMC_log_00.log",
            ),
            "r",
        ) as expected_file:
            expected_lines = expected_file.readlines()
            result_lines = result_file.readlines()
            for line_no, line in enumerate(expected_lines):
                floats_expected = []
                floats_result = []
                # Get just the floats
                for word in line.split():
                    try:
                        floats_expected.append(float(word))
                    except ValueError:
                        pass
                for word in result_lines[line_no].split():
                    try:
                        floats_result.append(float(word))
                    except ValueError:
                        pass
                # Skip the last float because it's just how long it took,
                # which will likely vary
                self.compare_equal(floats_expected[:-1], response=floats_result[:-1])

    def test_check_KMC_pickle_identical(self, run_simulation):
        output_morph_dir = run_simulation["output_parameter_dict"]["output_morph_dir"]
        morphology = run_simulation["output_parameter_dict"]["morphology"]
        KMC_dir = os.path.join(output_morph_dir, os.path.splitext(morphology)[0], "KMC")
        with open(
            os.path.join(
                TEST_ROOT,
                "assets",
                os.path.splitext(morphology)[0],
                "MKMC",
                "KMC",
                "KMC_results.pickle",
            ),
            "rb",
        ) as expected_pickle:
            expected = pickle.load(expected_pickle)
        with open(os.path.join(KMC_dir, "KMC_results.pickle"), "rb") as results_pickle:
            results = pickle.load(results_pickle)
        self.compare_equal(expected, response=results)


def teardown_module():
    shutil.rmtree(os.path.join(TEST_ROOT, "output_MKMC"))


if __name__ == "__main__":
    run_simulation()
