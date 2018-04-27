import os
import glob
import sys
import numpy as np
import random as R
import pickle
import subprocess as sp
from morphct.definitions import SINGLE_RUN_MOB_KMC_FILE
from morphct.code import helper_functions as hf


def main(AA_morphology_dict, CG_morphology_dict, CG_to_AAID_master, parameter_dict, chromophore_list):
    try:
        if parameter_dict['use_average_hop_rates']:
            print("Be advised: use_average_hop_rates is set to", repr(parameter_dict['use_average_hop_rates']) + ".")
            print("Orca-calculated energy levels will be ignored, and the following hop rates will be used:")
            print("Average Intra-molecular hop rate:", parameter_dict['average_intra_hop_rate'])
            print("Average Inter-molecular hop rate:", parameter_dict['average_inter_hop_rate'])
    except KeyError:
        pass
    # Attempt 2. PROGRAM THIS SERIALLY FIRST SO THAT IT WORKS
    # Determine the maximum simulation times based on the parameter dictionary
    simulation_times = parameter_dict['simulation_times']
    carrier_list = []
    # Modification: Rather than being clever here with the carriers, I'm just
    # going to create the master list of jobs that need running and then
    # randomly shuffle it. This will hopefully permit a similar number of holes
    # and electrons and lifetimes to be run simultaneously providing adequate
    # statistics more quickly
    for lifetime in simulation_times:
        for carrier_no in range(parameter_dict['number_of_holes_per_simulation_time']):
            carrier_list.append([carrier_no, lifetime, 'hole'])
        for carrier_no in range(parameter_dict['number_of_electrons_per_simulation_time']):
            carrier_list.append([carrier_no, lifetime, 'electron'])
    R.shuffle(carrier_list)
    proc_IDs = parameter_dict['proc_IDs']
    output_dir = parameter_dict['output_morph_dir'] + '/' + parameter_dict['morphology'][:-4] + '/KMC'
    jobs_list = [carrier_list[i:i + (int(np.ceil(len(carrier_list) / len(proc_IDs))))]
                 for i in range(0, len(carrier_list), int(np.ceil(len(carrier_list) / float(len(proc_IDs)))))]
    print("Writing job pickles for each CPU...")
    running_jobs = []
    for proc_ID, jobs in enumerate(jobs_list):
        pickle_name = output_dir + '/KMC_data_%02d.pickle' % (proc_ID)
        with open(pickle_name, 'wb+') as pickle_file:
            pickle.dump(jobs, pickle_file)
        print("KMC jobs for proc_ID", proc_ID, "written to KMC_data_%02d.pickle" % (proc_ID))
        # Open the required processes to execute the KMC jobs
        running_jobs.append(sp.Popen(['python', SINGLE_RUN_MOB_KMC_FILE, output_dir, str(proc_ID)]))
    # Wait for all jobs to complete
    [p.wait() for p in running_jobs]
    # Now combine all of the pickle files into one:
    print("All KMC jobs completed!")
    if parameter_dict['combine_KMC_results'] is True:
        print("Combining outputs...")
        combined_data = {}
        for proc_ID, jobs in enumerate(jobs_list):
            file_name = output_dir + '/KMC_results_%02d.pickle' % (proc_ID)
            # The pickle was repeatedly dumped to, in order to save time.
            # Each dump stream is self-contained, so iteratively unpickle to
            # add the new data.
            with open(file_name, 'rb') as pickle_file:
                pickled_data = pickle.load(pickle_file)
                for key, val in pickled_data.items():
                    if key not in combined_data:
                        combined_data[key] = val
                    else:
                        combined_data[key] += val
        # Write out the combined data
        with open(output_dir + '/KMC_results.pickle', 'wb+') as pickle_file:
            pickle.dump(combined_data, pickle_file)
        print("Complete data written to", output_dir + "/KMC_results.pickle.")
        print("Cleaning up...")
        # Delete any unneeded files
        for file_name in glob.glob(output_dir + '/KMC_results_*'):
            os.remove(file_name)
        for file_name in glob.glob(output_dir + '/KMC_slot_*'):
            os.remove(file_name)
    for file_name in glob.glob(output_dir + '/KMC_data*'):
        os.remove(file_name)
    return [AA_morphology_dict, CG_morphology_dict, CG_to_AAID_master, parameter_dict, chromophore_list]

if __name__ == '__main__':
    try:
        pickle_file = sys.argv[1]
    except:
        print("Please specify the pickle file to load to continue the pipeline from this point.")
    pickle_data = hf.load_pickle(pickle_file)
    AA_morphology_dict = pickle_data[0]
    CG_morphology_dict = pickle_data[1]
    CG_to_AAID_master = pickle_data[2]
    parameter_dict = pickle_data[3]
    chromophore_list = pickle_data[4]
    main(AA_morphology_dict, CG_morphology_dict, CG_to_AAID_master, parameter_dict, chromophore_list)
