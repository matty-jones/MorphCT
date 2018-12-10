import copy
import os
import pytest
import shutil
from morphct import run_MorphCT
from morphct.definitions import TEST_ROOT
from testing_tools import TestCommand
from morphct.code import helper_functions as hf


@pytest.fixture(
    scope="module",
    params=[
        {"voronoi": False, "hop_range": 10.0, "transiting": False, "koopmans": False},
        {"voronoi": True, "hop_range": 10.0, "transiting": False, "koopmans": False},
        {"voronoi": False, "hop_range": 5.0, "transiting": False, "koopmans": False},
        {"voronoi": False, "hop_range": 10.0, "transiting": True, "koopmans": False},
        {"voronoi": False, "hop_range": 10.0, "transiting": False, "koopmans": True},
    ],
)
def run_simulation(request):
    _voronoi = request.param["voronoi"]
    _hop_range_value = request.param["hop_range"]
    _transiting = request.param["transiting"]
    _koopmans = request.param["koopmans"]
    if _hop_range_value != 10.0:
        _hop_range = True
    else:
        _hop_range = False
    _modifiers = {"True": [], "False": []}
    for _key in request.param.keys():
        _boolean = eval("".join(["_", _key]))
        _modifiers[repr(_boolean)].append(_key)

    # ---==============================================---
    # ---======== Directory and File Structure ========---
    # ---==============================================---

    input_morph_dir = os.path.join(TEST_ROOT, "assets", "donor_polymer")
    output_morph_dir = os.path.join(TEST_ROOT, "output_OC")
    output_orca_dir = None
    input_device_dir = os.path.join(TEST_ROOT, "assets", "donor_polymer")
    output_device_dir = os.path.join(TEST_ROOT, "output_OC")

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
        True
    )  # Requires: Atomistic morphology, or molecular_dynamics
    execute_ZINDO = False  # Requires: obtain_chromophores
    execute_calculate_transfer_integrals = False  # Requires: execute_ZINDO
    execute_calculate_mobility = False  # Requires: calculate_transfer_integrals
    execute_device_simulation = (
        False
    )  # Requires: calculate_transfer_integrals for all device_components

    # ---==============================================---
    # ---============ Chromophore Parameters ==========---
    # ---==============================================---

    molecule_terminating_connections = {"C1": [[2, 1]], "C10": [[2, 1]]}
    AA_rigid_body_species = {}
    CG_site_species = {"A": "donor", "B": "none", "C": "none"}
    use_voronoi_neighbours = _voronoi
    maximum_hole_hop_distance = _hop_range
    maximum_electron_hop_distance = _hop_range
    permit_hops_through_opposing_chromophores = _transiting
    remove_orca_inputs = False
    remove_orca_outputs = True
    chromophore_length = 3

    # ---==============================================---
    # ---=== Chromophore Energy Scaling Parameters ====---
    # ---==============================================---

    chromophore_species = {
        "donor": {
            "literature_MO": -5.0,
            "target_DOS_std": 0.1,
            "reorganisation_energy": 0.3064,
            "species": "donor",
            "VRH_delocalisation": 2e-10,
        }
    }
    use_koopmans_approximation = _koopmans
    koopmans_hopping_prefactor = 1E-3

    # ---==============================================---
    # ---================= Begin run ==================---
    # ---==============================================---

    parameter_file = os.path.realpath(__file__)
    proc_IDs = hf.get_CPU_cores()
    parameter_names = [
        i
        for i in dir()
        if (not i.startswith("_"))
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
                "request",
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
            "OC",
            morphology.replace(".xml", "_post_run_HOOMD.pickle"),
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
    # Load the correct expected pickle
    pickle_name = os.path.join(
        input_morph_dir, "OC", morphology.replace(".xml", "_post_obtain_chromophores")
    )
    pickle_name_modified = "_".join([pickle_name] + sorted(_modifiers["True"]))
    expected_pickle_data = hf.load_pickle("".join([pickle_name_modified, ".pickle"]))
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
        for index, chromophore in enumerate(
            run_simulation["expected_chromophore_list"]
        ):
            for key in chromophore.__dict__.keys():
                self.compare_equal(
                    chromophore.__dict__[key],
                    response=run_simulation["output_chromophore_list"][index].__dict__[
                        key
                    ],
                )


def teardown_module():
    shutil.rmtree(os.path.join(TEST_ROOT, "output_OC"))


if __name__ == "__main__":

    class parameters:
        def __init__(self, param):
            self.param = param

    run_simulation(
        parameters(
            {"voronoi": False, "hop_range": 10.0, "transiting": False, "koopmans": True}
        )
    )
