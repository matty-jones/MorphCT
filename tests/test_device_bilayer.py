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

    input_morph_dir = os.path.join(TEST_ROOT, "output_DBL", "morph_inputs")
    output_morph_dir = os.path.join(TEST_ROOT, "output_DBL", "morph_outputs")
    output_orca_dir = None
    input_device_dir = os.path.join(TEST_ROOT, "output_DBL", "device_inputs")
    output_device_dir = os.path.join(TEST_ROOT, "output_DBL", "device_outputs")

    # ---==============================================---
    # ---========== Input Morphology Details ==========---
    # ---==============================================---

    morphology = None
    input_sigma = 1.0
    device_morphology = "bilayer"
    device_components = {
        0: "donor_crystal",
        1: "acceptor_crystal",
        2: "mixed_crystal_bilayer",
    }
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
    execute_calculate_mobility = False  # Requires: calculate_transfer_integrals
    execute_device_simulation = (
        True
    )  # Requires: calculate_transfer_integrals for all device_components

    # ---==============================================---
    # ---=== General Kinetic Monte Carlo Parameters ===---
    # ---==============================================---

    # ---=== Universal KMC Parameters ===---
    system_temperature = 290
    use_simple_energetic_penalty = False
    record_carrier_history = True
    use_VRH = True

    # ---=== Device Specific KMC Parameters ===---
    # Material parameters
    absorption_coefficient = 1.3E4
    relative_permittivity = 3
    donor_HOMO = -5.3
    acceptor_LUMO = -3.9
    cathode_work_function = -4.2
    anode_work_function = -5.0
    recombination_rate = 1E9
    coulomb_capture_radius = 1E-9
    wrap_device_xy = False
    exciton_lifetime = 0.5E-9
    forster_radius = 4.3E-9
    hopping_prefactor = 1E-4
    MA_prefactor = 1E11
    MA_localisation_radius = 1E-9

    # External parameters
    incident_flux = 100
    incident_wavelength = 500E-9

    # Simulation
    voltage_sweep = [-0.1]
    morphology_cell_size = 3E-8
    minimum_number_of_photoinjections = 10
    fastest_event_allowed = 1E-18
    slowest_event_allowed = 1E-8
    disable_dark_injection = True
    disable_coulombic = True
    output_log_to_stdout = True

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
        shutil.rmtree(os.path.join(TEST_ROOT, "output_DBL"))
    except OSError:
        pass

    shutil.copytree(
        os.path.join(TEST_ROOT, "assets", "device_bilayer"),
        os.path.join(TEST_ROOT, "output_DBL"),
    )

    run_MorphCT.simulation(
        **parameters
    )  # Execute MorphCT using these simulation parameters


# ---==============================================---
# ---================= Run Tests ==================---
# ---==============================================---


class TestCompareOutputs(TestCommand):
    def test_wobbey(self, run_simulation):
        return True


def teardown_module():
    shutil.rmtree(os.path.join(TEST_ROOT, "output_DBL"))


if __name__ == "__main__":
    run_simulation()
