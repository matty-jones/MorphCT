import distutils.spawn
import os
import pickle
import sys
import subprocess as sp
import time as T
from morphct.code import helper_functions as hf


if __name__ == '__main__':
    morphology_file = sys.argv[1]
    CPU_rank = int(sys.argv[2])
    overwrite = False
    try:
        overwrite = bool(sys.argv[3])
    except:
        pass
    morphology_name = morphology_file[hf.find_index(morphology_file, '/')[-1] + 1:]
    orca_path = distutils.spawn.find_executable("orca")
    input_dir = morphology_file + '/chromophores/input_orca'
    log_file = input_dir.replace('/input_orca', '/orca_log_' + str(CPU_rank) + '.log')
    output_dir = morphology_file + '/chromophores/output_orca'
    pickle_file_name = input_dir.replace('input_orca', 'orca_jobs.pickle')
    with open(pickle_file_name, 'rb') as pickle_file:
        jobs_list = pickle.load(pickle_file)
    jobs_to_run = jobs_list[CPU_rank]
    hf.write_to_file(log_file, ["Found " + str(len(jobs_to_run)) + " jobs to run."])
    t0 = T.time()
    for job in jobs_to_run:
        t1 = T.time()
        hf.write_to_file(log_file, ["Running job " + str(job) + "..."])
        output_file_name = job.replace('.inp', '.out').replace('input_orca', 'output_orca')
        # Check if file exists already
        if overwrite is False:
            try:
                with open(output_file_name, 'r') as test_file:
                    pass
                hf.write_to_file(log_file, [output_file_name + " already exists, and"
                                            " overwrite_current_data is " + repr(overwrite)
                                            + "! skipping..."])
                continue
            except IOError:
                pass
        orca_job = sp.Popen([str(orca_path), str(job)], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
        job_PID = orca_job.pid
        # Taskset stuff
        # try:
        #     affinity_job = sp.Popen(['taskset', '-pc', str(CPU_rank), str(job_PID)],
        #                             stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
        #     # hf.write_to_file(log_file, affinity_job[0].split('\n'))  # stdOut for affinity set
        #     # hf.write_to_file(log_file, affinity_job[1].split('\n'))  # stdErr for affinity set
        # except OSError:
        #     hf.write_to_file(log_file, ["Taskset command not found, skipping setting of processor affinities..."])
        orca_shell_output = orca_job.communicate()
        # Write the outputFile:
        hf.write_to_file(output_file_name, orca_shell_output[0].decode().split('\n'), mode='output_file')
        hf.write_to_file(log_file, orca_shell_output[1].decode().split('\n'))  # std_err
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
        elapsed_time = '%.1f' % (float(elapsed_time))
        hf.write_to_file(log_file, ["Job " + str(job) + " completed in " + elapsed_time + " " + time_units])
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
    elapsed_time = '%.1f' % (float(elapsed_time))
    hf.write_to_file(log_file, ["All jobs completed in " + elapsed_time + " " + time_units])
    hf.write_to_file(log_file, ["Exiting normally..."])
