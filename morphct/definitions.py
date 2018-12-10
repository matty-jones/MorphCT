import os

PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
TEST_ROOT = os.path.join(os.path.dirname(PROJECT_ROOT), "tests")
SINGLE_ORCA_RUN_FILE = os.path.join(PROJECT_ROOT, "code", "single_core_run_orca.py")
SINGLE_RUN_MOB_KMC_FILE = os.path.join(
    PROJECT_ROOT, "code", "single_core_run_mob_KMC.py"
)
SINGLE_RUN_DEVICE_KMC_FILE = os.path.join(
    PROJECT_ROOT, "code", "single_core_run_device_KMC.py"
)
