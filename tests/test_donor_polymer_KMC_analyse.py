import csv
import os
import pytest
import shutil
import subprocess
from morphct import run_MorphCT
from morphct.definitions import TEST_ROOT
from testing_tools import TestCommand
from morphct.code import helper_functions as hf


@pytest.fixture(
    scope="module",
    params=[
        "",
        "-t",
        '-x "TEST_OUTPUT"',
        "-sd 0.0,1.0",
        "COMBINE_KMC",
        "-crd auto -cod auto -cfd auto",
        "-ctd 0.1 -crd 1.0 -cod 10.0",
    ],
)
def run_simulation(request):
    flags = request.param

    # ---==============================================---
    # ---=============== Setup Prereqs ================---
    # ---==============================================---

    output_dir = os.path.join(TEST_ROOT, "output_KMCA")

    try:
        shutil.rmtree(output_dir)
    except OSError:
        pass
    os.makedirs(os.path.join(output_dir, "code"))
    os.makedirs(os.path.join(output_dir, "figures"))
    shutil.copy(
        os.path.join(
            TEST_ROOT,
            "assets",
            "donor_polymer",
            "MKMC",
            "donor_polymer_post_calculate_mobility.pickle",
        ),
        os.path.join(
            output_dir, "code", "".join([os.path.split(output_dir)[1], ".pickle"])
        ),
    )
    if flags != "COMBINE_KMC":
        shutil.copytree(
            os.path.join(TEST_ROOT, "assets", "donor_polymer", "MKMC", "KMC"),
            os.path.join(output_dir, "KMC"),
        )
        command = [
            component
            for component in ["KMCAnalyse"]
            + flags.split()
            + ["-b", "AGG"]
            + [output_dir]
        ]
        if "-s" in flags:
            command.append(output_dir)
    else:
        shutil.copytree(
            os.path.join(TEST_ROOT, "assets", "donor_polymer", "MKMC", "KMC2"),
            os.path.join(output_dir, "KMC"),
        )
        command = [component for component in ["KMCAnalyse", "-b AGG", output_dir]]
    print("Executing command", command)
    subprocess.Popen(command).communicate()
    return flags


# ---==============================================---
# ---================= Run Tests ==================---
# ---==============================================---


class TestCompareOutputs(TestCommand):
    def test_check_results_exists(self, run_simulation):
        self.confirm_file_exists(os.path.join(TEST_ROOT, "output_KMCA", "results.csv"))

    def test_compare_results(self, run_simulation):
        flags = ""
        if run_simulation == "COMBINE_KMC":
            return True
        elif len(run_simulation) > 0:
            flags = run_simulation.split()[0]
        with open(
            os.path.join(
                TEST_ROOT,
                "assets",
                "donor_polymer",
                "KMCA",
                "".join(["results", flags, ".csv"]),
            ),
            "r",
        ) as expected_csv:
            reader = csv.reader(expected_csv)
            expected_dict = {rows[0]: rows[1] for rows in reader}
        with open(
            os.path.join(TEST_ROOT, "output_KMCA", "results.csv"), "r"
        ) as results_csv:
            reader = csv.reader(results_csv)
            results_dict = {rows[0]: rows[1] for rows in reader}
        for key in expected_dict:
            try:
                expected_dict[key] = float(expected_dict[key])
            except ValueError:
                pass
        for key in results_dict:
            try:
                results_dict[key] = float(results_dict[key])
            except ValueError:
                pass
        self.compare_equal(expected_dict, response=results_dict)

    def test_check_network_figure(self, run_simulation):
        if "-t" in run_simulation:
            self.confirm_file_exists(
                os.path.join(TEST_ROOT, "output_KMCA", "figures", "01_3d_hole.png")
            )
        else:
            self.confirm_file_exists(
                os.path.join(TEST_ROOT, "output_KMCA", "figures", "01_3d_hole.png"),
                negate=True,
            )

    def test_check_stack_figure(self, run_simulation):
        if "-t" in run_simulation:
            self.confirm_file_exists(
                os.path.join(TEST_ROOT, "output_KMCA", "figures", "03_clusters.png")
            )
        else:
            self.confirm_file_exists(
                os.path.join(TEST_ROOT, "output_KMCA", "figures", "03_clusters.png"),
                negate=True,
            )

    def test_check_neighbour_hist_figure(self, run_simulation):
        self.confirm_file_exists(
            os.path.join(
                TEST_ROOT, "output_KMCA", "figures", "04_neighbour_hist_donor.png"
            )
        )

    def test_check_delta_E_figure(self, run_simulation):
        self.confirm_file_exists(
            os.path.join(TEST_ROOT, "output_KMCA", "figures", "06_donor_delta_E_ij.png")
        )

    def test_check_anisotropy_figure(self, run_simulation):
        if "-t" in run_simulation:
            self.confirm_file_exists(
                os.path.join(
                    TEST_ROOT, "output_KMCA", "figures", "08_anisotropy_hole.png"
                )
            )
        else:
            self.confirm_file_exists(
                os.path.join(
                    TEST_ROOT, "output_KMCA", "figures", "08_anisotropy_hole.png"
                ),
                negate=True,
            )

    def test_check_transfer_integral_mols_figure(self, run_simulation):
        self.confirm_file_exists(
            os.path.join(
                TEST_ROOT,
                "output_KMCA",
                "figures",
                "10_donor_transfer_integral_mols.png",
            )
        )

    def test_check_transfer_integral_stacks_figure(self, run_simulation):
        self.confirm_file_exists(
            os.path.join(
                TEST_ROOT,
                "output_KMCA",
                "figures",
                "12_donor_transfer_integral_clusters.png",
            )
        )

    def test_check_hopping_rate_mols_figure(self, run_simulation):
        self.confirm_file_exists(
            os.path.join(
                TEST_ROOT, "output_KMCA", "figures", "14_donor_hopping_rate_mols.png"
            )
        )

    def test_check_hopping_rate_stacks_figure(self, run_simulation):
        self.confirm_file_exists(
            os.path.join(
                TEST_ROOT,
                "output_KMCA",
                "figures",
                "16_donor_hopping_rate_clusters.png",
            )
        )

    def test_check_lin_MSD_figure(self, run_simulation):
        self.confirm_file_exists(
            os.path.join(TEST_ROOT, "output_KMCA", "figures", "18_lin_MSD_hole.png")
        )

    def test_check_semi_log_MSD_figure(self, run_simulation):
        self.confirm_file_exists(
            os.path.join(
                TEST_ROOT, "output_KMCA", "figures", "20_semi_log_MSD_hole.png"
            )
        )

    def test_check_log_MSD_figure(self, run_simulation):
        self.confirm_file_exists(
            os.path.join(TEST_ROOT, "output_KMCA", "figures", "22_log_MSD_hole.png")
        )

    def test_check_log_MSD_figure(self, run_simulation):
        self.confirm_file_exists(
            os.path.join(
                TEST_ROOT, "output_KMCA", "figures", "24_total_hop_freq_hole.png"
            )
        )

    def test_check_log_MSD_figure(self, run_simulation):
        self.confirm_file_exists(
            os.path.join(
                TEST_ROOT, "output_KMCA", "figures", "26_net_hop_freq_hole.png"
            )
        )

    def test_check_log_MSD_figure(self, run_simulation):
        self.confirm_file_exists(
            os.path.join(
                TEST_ROOT, "output_KMCA", "figures", "28_hop_discrepancy_hole.png"
            )
        )

    def test_check_log_MSD_figure(self, run_simulation):
        self.confirm_file_exists(
            os.path.join(
                TEST_ROOT, "output_KMCA", "figures", "30_hole_displacement_dist.png"
            )
        )

    def test_check_log_MSD_figure(self, run_simulation):
        self.confirm_file_exists(
            os.path.join(
                TEST_ROOT, "output_KMCA", "figures", "32_hole_cluster_dist.png"
            )
        )

    def test_check_log_MSD_figure(self, run_simulation):
        self.confirm_file_exists(
            os.path.join(
                TEST_ROOT, "output_KMCA", "figures", "34_orientation_hist_donor.png"
            )
        )

    def test_check_hop_vec_figure(self, run_simulation):
        if "-t" in run_simulation:
            self.confirm_file_exists(
                os.path.join(TEST_ROOT, "output_KMCA", "figures", "36_hop_vec_hole.png")
            )
        else:
            self.confirm_file_exists(
                os.path.join(
                    TEST_ROOT, "output_KMCA", "figures", "36_hop_vec_hole.png"
                ),
                negate=True,
            )

    def test_check_anisotropy_sequence_figure(self, run_simulation):
        if "-s" in run_simulation:
            self.confirm_file_exists(os.path.join(os.getcwd(), "anisotropy_hole.png"))

    def test_check_mobility_sequence_figure(self, run_simulation):
        if "-s" in run_simulation:
            self.confirm_file_exists(os.path.join(os.getcwd(), "mobility_hole.png"))


def teardown_module():
    shutil.rmtree(os.path.join(TEST_ROOT, "output_KMCA"))
    os.remove(os.path.join(os.getcwd(), "anisotropy_hole.png"))
    os.remove(os.path.join(os.getcwd(), "mobility_hole.png"))


if __name__ == "__main__":

    class parameters:
        def __init__(self, param):
            self.param = param

    run_simulation(parameters("-t"))
