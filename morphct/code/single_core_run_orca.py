import distutils.spawn
import os
import pickle
import sys
import subprocess as sp
import time as T
from morphct.code import helper_functions as hf


### There are still so many edge cases that the following function doesn't always
### work. Use the one below it instead.
# def check_job_output(orca_lines):
#     # Sometimes ORCA terminates normally even though the SCF didn't converge.
#     # Also, sometimes ORCA will fail at the SCF but not use the word "ERROR"
#     for line in orca_lines:
#         # Any "ERROR" lines should come first
#         if "ERROR" in line:
#             return False
#         # Only if there are no "ERROR" lines and "ORCA TERMINATED NORMALLY" is
#         # near/at the bottom was the job successful.
#         if "ORCA TERMINATED NORMALLY" in line:
#             return True
#     # If no "ERROR" but also no "ORCA TERMINATED NORMALLY", then the job failed.
#     return False


def check_job_output(data_file, job):
    record_MO_data = False
    orbital_data = []
    for line in data_file:
        if "ORBITAL ENERGIES" in line:
            # Next line begins the MO data
            record_MO_data = True
            continue
        if record_MO_data is True:
            if "MOLECULAR ORBITALS" in line:
                # Don't need anything else from the output file
                break
            data_in_line = []
            for element in line.split(" "):
                if len(element) > 1:
                    try:
                        data_in_line.append(float(element))
                    except ValueError:
                        continue
            if len(data_in_line) == 4:
                orbital_data.append(data_in_line)
    for i in range(len(orbital_data)):
        if orbital_data[i][1] == 0:
            # This line is the first unoccupied orbital - i.e. LUMO
            LUMO = orbital_data[i][3]
            HOMO = orbital_data[i - 1][3]
            HOMO_1 = orbital_data[i - 2][3]
            LUMO_1 = orbital_data[i + 1][3]
            # Don't need any other orbitals
            break
    if record_MO_data is False:
        # Molecular orbital data not present in this file
        print("Molecular orbital data not present for", job, ". Skipping deletion")
        return False
    print("Molecular orbital data detected. Deleting", job, "...")
    return True


def main():
    morph_orca_dir = sys.argv[1]
    morph_output_dir = sys.argv[2]
    CPU_rank = int(sys.argv[3])
    overwrite = bool(int(sys.argv[4]))
    delete_inputs = bool(int(sys.argv[5]))
    morphology_name = os.path.split(morph_orca_dir[1])
    try:
        orca_path = os.environ["ORCA_BIN"]
    except KeyError:
        orca_path = distutils.spawn.find_executable("orca")
    input_dir = os.path.join(morph_orca_dir, "chromophores", "input_orca")
    log_file = os.path.join(
        morph_output_dir, "chromophores", "orca_log_{:02d}.log".format(CPU_rank)
    )
    output_dir = os.path.join(morph_orca_dir, "chromophores", "output_orca")
    pickle_file_name = input_dir.replace("input_orca", "orca_jobs.pickle")
    with open(pickle_file_name, "rb") as pickle_file:
        jobs_list = pickle.load(pickle_file)
    jobs_to_run = jobs_list[CPU_rank]
    hf.write_to_file(log_file, ["Found {:d} jobs to run.".format(len(jobs_to_run))])
    t0 = T.time()
    for job in jobs_to_run:
        t1 = T.time()
        hf.write_to_file(log_file, ["Running job {:s}...".format(job)])
        output_file_name = job.replace(".inp", ".out").replace(
            "input_orca", "output_orca"
        )
        # Check if file exists already
        if overwrite is False:
            try:
                with open(output_file_name, "r") as test_file:
                    pass
                hf.write_to_file(
                    log_file,
                    [
                        "".join(
                            [
                                output_file_name,
                                " already exists, and overwrite_current_data is ",
                                repr(overwrite),
                                "! skipping...",
                            ]
                        )
                    ],
                )
                continue
            except IOError:
                pass
        orca_job = sp.Popen(
            [str(orca_path), str(job)], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE
        )
        job_PID = orca_job.pid
        # Taskset stuff
        # try:
        #     affinity_job = sp.Popen(['taskset', '-pc', str(CPU_rank), str(job_PID)],
        #                             stdin=sp.PIPE, stdout=sp.PIPE,
        #                             stderr=sp.PIPE).communicate()
        #     # hf.write_to_file(log_file, affinity_job[0].split('\n'))
        #     # hf.write_to_file(log_file, affinity_job[1].split('\n'))
        # except OSError:
        #     hf.write_to_file(log_file, ["Taskset command not found, skipping setting"
        #                                 " of processor affinities..."])
        orca_shell_output = orca_job.communicate()
        orca_stdout = orca_shell_output[0].decode().split("\n")
        orca_stderr = orca_shell_output[1].decode().split("\n")
        # Write the outputFile:
        try:
            hf.write_to_file(output_file_name, orca_stdout, mode="output_file")
        except OSError as e1:
            file_written = False
            print(e1)
            for attempt in range(5):
                print("Waiting for 10 minutes to try and write again...")
                T.sleep(600)
                try:
                    hf.write_to_file(output_file_name, orca_stdout, mode="output_file")
                    file_written = True
                except OSError as e2:
                    print(e2)
                    continue
            if not file_written:
                print("--== CRITICAL ERROR ==--")
                print("Unable to write file due to above errors. Please ammend"
                      " appropriately and resubmit your job")
                if parameter_dict["remove_orca_inputs"]:
                    print("Deleting the remaining orca inputs...")
                    shutil.rmtree(orca_input_dir)
                if parameter_dict["remove_orca_outputs"]:
                    print("Deleting the remaining orca inputs...")
                    shutil.rmtree(orca_output_dir)
                print("Exiting...")
                exit()
        if delete_inputs:
            output_ok = check_job_output(orca_stdout, job)
            if output_ok:
                hf.write_to_file(
                    log_file,
                    [
                        "Output OK and remove_orca_inputs set.",
                        "".join(["Deleting ", os.path.splitext(job)[0], " inputs..."]),
                    ],
                )
                for extension in [
                    ".inp",
                    ".gbw",
                    ".prop",
                    ".tmp",
                    ".ges",
                    "_property.txt",
                ]:
                    try:
                        os.remove(job.replace(".inp", extension))
                    except FileNotFoundError:
                        # Already deleted
                        pass
        if len(orca_stderr) > 1:
            # Write any errors
            hf.write_to_file(log_file, orca_stderr)
        t2 = T.time()
        elapsed_time = float(t2) - float(t1)
        if elapsed_time < 60:
            time_units = "seconds."
        elif elapsed_time < 3600:
            elapsed_time /= 60.0
            time_units = "minutes."
        elif elapsed_time < 86400:
            elapsed_time /= 3600.0
            time_units = "hours."
        else:
            elapsed_time /= 86400.0
            time_units = "days."
        hf.write_to_file(
            log_file,
            [
                "Job {0:s} completed in {1:.2f} {2:s}\n".format(
                    job, elapsed_time, time_units
                )
            ],
        )
        # Now check the output file and delete the input files if we don't need them
        # any more
    t3 = T.time()
    elapsed_time = float(t3) - float(t0)
    if elapsed_time < 60:
        time_units = "seconds."
    elif elapsed_time < 3600:
        elapsed_time /= 60.0
        time_units = "minutes."
    elif elapsed_time < 86400:
        elapsed_time /= 3600.0
        time_units = "hours."
    else:
        elapsed_time /= 86400.0
        time_units = "days."
    hf.write_to_file(
        log_file,
        ["All jobs completed in {0:.2f} {1:s}".format(elapsed_time, time_units)],
    )
    hf.write_to_file(log_file, ["Exiting normally..."])


if __name__ == "__main__":
    main()
