import copy
import os
import pytest
import shutil
import sys
from morphct import run_MorphCT
from morphct.definitions import TEST_ROOT
from testing_tools import TestCommand
from morphct.code import helper_functions as hf
from morphct.code.helper_functions import has_hoomd


@pytest.fixture(scope="module")
def run_simulation():
    # ---==============================================---
    # ---======== Directory and File Structure ========---
    # ---==============================================---

    input_morph_dir = os.path.join(TEST_ROOT, "assets", "donor_polymer")
    output_morph_dir = os.path.join(TEST_ROOT, "output_RH")
    output_orca_dir = None
    input_device_dir = os.path.join(TEST_ROOT, "assets", "donor_polymer")
    output_device_dir = os.path.join(TEST_ROOT, "output_RH")

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
    execute_molecular_dynamics = True  # Requires: fine_graining
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
    # ---=========== Forcefield Parameters ============---
    # ---==============================================---

    CG_to_template_dirs = {
        "A": os.path.join(TEST_ROOT, "assets", "donor_polymer"),
        "B": os.path.join(TEST_ROOT, "assets", "donor_polymer"),
        "C": os.path.join(TEST_ROOT, "assets", "donor_polymer"),
    }
    CG_to_template_force_fields = {
        "A": "test_FF.xml",
        "B": "test_FF.xml",
        "C": "test_FF.xml",
    }
    pair_r_cut = 10.0
    pair_dpd_gamma_val = 0.0

    # ---==============================================---
    # ---===== Molecular Dynamics Phase Parameters ====---
    # ---==============================================---

    number_of_phases = 4
    temperatures = [1.0]
    taus = [1.0]
    pair_types = ["none", "dpd", "lj", "lj", "lj", "lj", "lj", "lj"]
    bond_types = ["harmonic"]
    angle_types = ["harmonic"]
    dihedral_types = ["opls"]
    integration_targets = ["all"]
    timesteps = [1E-3, 1E-3, 1E-7, 1E-4]
    durations = [1E5, 5E4, 1E4, 5E4]
    termination_conditions = ["ke_min", "max_t", "max_t", "max_t"]
    group_anchorings = ["all", "all", "all", "none"]
    dcd_file_write = True
    dcd_file_dumpsteps = [0]

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
            "RH",
            morphology.replace(".xml", "_post_fine_graining.pickle"),
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
            input_morph_dir, "RH", morphology.replace(".xml", "_post_run_HOOMD.pickle")
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


@pytest.mark.skipif(
    not has_hoomd, reason="HOOMD 1.3 is not installed on this system."
)
class TestCompareOutputs(TestCommand):
    def test_check_AA_morphology_dict_len(self, run_simulation):
        for key in run_simulation["expected_AA_morphology_dict"]:
            try:
                self.compare_equal(
                    len(run_simulation["expected_AA_morphology_dict"][key]),
                    response=len(run_simulation["output_AA_morphology_dict"][key]),
                )
            except TypeError:
                self.compare_equal(
                    run_simulation["expected_AA_morphology_dict"][key],
                    response=run_simulation["output_AA_morphology_dict"][key],
                )

    @pytest.mark.skipif(
        sys.platform != "darwin",
        reason="Expected output will change on non-OSX platforms.",
    )
    def test_check_AA_morphology_dict_direct(self, run_simulation):
        for key in run_simulation["expected_AA_morphology_dict"]:
            self.compare_equal(
                run_simulation["expected_AA_morphology_dict"][key],
                response=run_simulation["output_AA_morphology_dict"][key],
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

    def test_check_morphology_phase0_created(self, run_simulation):
        output_morph_dir = run_simulation["output_parameter_dict"]["output_morph_dir"]
        morphology = run_simulation["output_parameter_dict"]["morphology"]
        morph_dir = os.path.join(
            output_morph_dir, os.path.splitext(morphology)[0], "morphology"
        )
        self.confirm_file_exists(os.path.join(morph_dir, "phase0_" + morphology))

    def test_check_morphology_phase1_created(self, run_simulation):
        output_morph_dir = run_simulation["output_parameter_dict"]["output_morph_dir"]
        morphology = run_simulation["output_parameter_dict"]["morphology"]
        morph_dir = os.path.join(
            output_morph_dir, os.path.splitext(morphology)[0], "morphology"
        )
        self.confirm_file_exists(os.path.join(morph_dir, "phase1_" + morphology))
        self.confirm_file_exists(
            os.path.join(morph_dir, "phase1_" + morphology.replace("xml", "dcd"))
        )

    def test_check_morphology_phase2_created(self, run_simulation):
        output_morph_dir = run_simulation["output_parameter_dict"]["output_morph_dir"]
        morphology = run_simulation["output_parameter_dict"]["morphology"]
        morph_dir = os.path.join(
            output_morph_dir, os.path.splitext(morphology)[0], "morphology"
        )
        self.confirm_file_exists(os.path.join(morph_dir, "phase2_" + morphology))
        self.confirm_file_exists(
            os.path.join(morph_dir, "phase2_" + morphology.replace("xml", "dcd"))
        )

    def test_check_morphology_phase3_created(self, run_simulation):
        output_morph_dir = run_simulation["output_parameter_dict"]["output_morph_dir"]
        morphology = run_simulation["output_parameter_dict"]["morphology"]
        morph_dir = os.path.join(
            output_morph_dir, os.path.splitext(morphology)[0], "morphology"
        )
        self.confirm_file_exists(os.path.join(morph_dir, "phase3_" + morphology))
        self.confirm_file_exists(
            os.path.join(morph_dir, "phase3_" + morphology.replace("xml", "dcd"))
        )

    def test_check_morphology_phase4_created(self, run_simulation):
        output_morph_dir = run_simulation["output_parameter_dict"]["output_morph_dir"]
        morphology = run_simulation["output_parameter_dict"]["morphology"]
        morph_dir = os.path.join(
            output_morph_dir, os.path.splitext(morphology)[0], "morphology"
        )
        self.confirm_file_exists(os.path.join(morph_dir, "phase4_" + morphology))
        self.confirm_file_exists(
            os.path.join(morph_dir, "phase4_" + morphology.replace("xml", "dcd"))
        )

    def test_check_morphology_final_created(self, run_simulation):
        output_morph_dir = run_simulation["output_parameter_dict"]["output_morph_dir"]
        morphology = run_simulation["output_parameter_dict"]["morphology"]
        morph_dir = os.path.join(
            output_morph_dir, os.path.splitext(morphology)[0], "morphology"
        )
        self.confirm_file_exists(os.path.join(morph_dir, "final_" + morphology))
        self.confirm_file_exists(
            os.path.join(morph_dir, "energies_" + morphology.replace("xml", "log"))
        )

    def test_check_morphology_final_output_pickle(self, run_simulation):
        output_morph_dir = run_simulation["output_parameter_dict"]["output_morph_dir"]
        morphology = run_simulation["output_parameter_dict"]["morphology"]
        morph_dir = os.path.join(
            output_morph_dir, os.path.splitext(morphology)[0], "morphology"
        )
        final_MD_output = hf.load_morphology_xml(
            os.path.join(morph_dir, "final_" + morphology)
        )
        self.compare_equal(
            final_MD_output, response=run_simulation["output_AA_morphology_dict"]
        )


def teardown_module():
    shutil.rmtree(os.path.join(TEST_ROOT, "output_RH"))


if __name__ == "__main__":
    run_simulation()
