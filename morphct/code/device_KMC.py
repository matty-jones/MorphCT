import os
import pickle
import sys
import numpy as np
import subprocess as sp
from morphct.definitions import SINGLE_RUN_DEVICE_KMC_FILE
from morphct.code import helper_functions as hf


class morphology_moiety:
    def __init__(self, mol_morph_name, parameter_dict):
        chromophore_list_location = os.path.join(
            parameter_dict["output_morph_dir"],
            mol_morph_name,
            "code",
            "".join([mol_morph_name, ".pickle"]),
        )
        pickle_data = hf.load_pickle(chromophore_list_location)
        self.AA_morphology_dict = pickle_data[0]
        self.parameter_dict = pickle_data[3]
        self.chromophore_list = pickle_data[4]
        self.carrier_type = self.get_carrier_type()
        # Now add the occupation data to the chromophoreLists so that we can
        # prevent double occupation in the simulations.
        # The occupied property is a list that contains the device moiety
        # coordinates where the chromophore is occupied.
        for index, chromophore in enumerate(self.chromophore_list):
            chromophore.occupied = []

    def get_carrier_type(self):
        species_present = []
        for chromophore in self.chromophore_list:
            species_present.append(chromophore.species)
        if len(set(species_present)) == 1:
            if species_present[0].lower() == "donor":
                return "hole"
            elif species_present[0].lower() == "acceptor":
                return "electron"
            else:
                print("Error in chromophore:")
                for key, val in chromophore.__dict__.items():
                    print(key, "=", val)
                raise SystemError("Chromophore species is neither donor nor acceptor")
        else:
            return "both"


class chromophore_data_container:
    # A helper class that contains all of the chromophore data for ease of
    # access from anywhere
    def __init__(self, device_array, moiety_dictionary, wrapxy):
        self.device_array = device_array
        self.moiety_dictionary = moiety_dictionary
        self.wrapxy = wrapxy

    def return_chromophore_list(self, device_position):
        device_moiety_type = self.device_array[tuple(device_position)]
        return self.moiety_dictionary[device_moiety_type].chromophore_list

    def return_specific_chromophore(self, device_position, chromo_ID):
        device_moiety_type = self.device_array[tuple(device_position)]
        return self.moiety_dictionary[device_moiety_type].chromophore_list[chromo_ID]

    def return_random_chromophore(self, device_position):
        device_moiety_type = self.device_array[tuple(device_position)]
        return np.random.choice(
            self.moiety_dictionary[device_moiety_type].chromophore_list
        )

    def return_closest_chromophore_to_position(self, device_position, desired_position):
        closest_chromo_ID = None
        # Check that there is an eligible device position that exists at these
        # coordinates (i.e. there is a hop taking place within the active layer)
        # Find out which axis is out of index
        for axis_no, val in enumerate(device_position):
            if val >= self.device_array.shape[axis_no]:
                if axis_no == 2:
                    # Leaving out of the top of the device
                    return "top"
                elif self.wrapxy:
                    # Bring it in on the reverse side
                    device_position[axis_no] = 0
                else:
                    return "out of bounds"
            if val < 0:
                if axis_no == 2:
                    # Leaving out of the bottom of the device
                    return "bottom"
                elif self.wrapxy:
                    # Bring it in on the reverse side
                    device_position[axis_no] = self.device_array.shape[axis_no] - 1
                else:
                    return "out of bounds"
        device_moiety_type = self.device_array[tuple(device_position)]
        positions = np.array(
            [
                chromo.posn
                for chromo in self.moiety_dictionary[
                    device_moiety_type
                ].chromophore_list
            ]
        )
        distances = np.sqrt(
            np.sum((positions - np.array(desired_position)) ** 2, axis=1)
        )
        closest_chromo_ID = np.argmin(distances)
        return self.moiety_dictionary[device_moiety_type].chromophore_list[
            closest_chromo_ID
        ]


class morphology_data_container:
    # A helper class that contains all of the chromophore data for ease of
    # access from anywhere
    def __init__(self, device_array, moiety_dictionary):
        self.device_array = device_array
        self.moiety_dictionary = moiety_dictionary

    def return_AA_morphology(self, device_position):
        device_moiety_type = self.device_array[tuple(device_position)]
        return self.moiety_dictionary[device_moiety_type].AA_morphology_dict

    def return_device_moiety_type(self, device_position):
        return self.device_array[tuple(device_position)]


def load_device_morphology(parameter_dict):
    device_dir = os.path.join(
        parameter_dict["input_device_dir"], parameter_dict["device_morphology"]
    )
    y_slices = os.listdir(device_dir)
    # Initialize the array of the correct size (assumes cubic morphology)
    device_array = np.zeros([len(y_slices)] * 3, dtype=int)
    for y_val, file_name in enumerate(y_slices):
        # Load the ySlice as-presented in the input files
        y_slice = np.loadtxt(os.path.join(device_dir, file_name), dtype=int)
        if len(y_slice.shape) > 0:
            # The z-origin is at the top, and we need it at the bottom, so turn
            # the array upside down
            y_slice = np.flipud(y_slice)
            # Now populate the array
            for z_val, z_row in enumerate(y_slice):
                for x_val, datum in enumerate(z_row):
                    device_array[x_val, y_val, z_val] = datum
        else:
            # Can't flipud and iterate over a zero-length array (one number), so
            # assign it this way instead.
            device_array[0, y_val, 0] = int(y_slice)
    moiety_dictionary = {}
    for moiety_ID in np.unique(device_array):
        moiety_dictionary[moiety_ID] = morphology_moiety(
            parameter_dict["device_components"][moiety_ID], parameter_dict
        )
    return device_array, moiety_dictionary


def main(parameter_dict):
    # Get the random seed now for all the child processes
    if parameter_dict["random_seed_override"] is not None:
        np.random.seed(parameter_dict["random_seed_override"])
    # First job will be to load in the device morphology, when I work out what
    # format I want it to be.
    device_array, moiety_dictionary = load_device_morphology(parameter_dict)
    # Initialise the helperClass to obtain all of the chromophoreData required,
    # allowing it be accessed globally
    chromophore_data = chromophore_data_container(
        device_array, moiety_dictionary, parameter_dict["wrap_device_xy"]
    )
    morphology_data = morphology_data_container(device_array, moiety_dictionary)
    # Write these classes out to a pickle file so that they can be loaded by the
    # child processes later
    to_pickle = [device_array, chromophore_data, morphology_data, parameter_dict]
    save_directory = os.path.join(
        parameter_dict["output_device_dir"], parameter_dict["device_morphology"], "code"
    )
    if parameter_dict["overwrite_current_data"] is True:
        with open(
            os.path.join(save_directory, "device_data.pickle"), "wb+"
        ) as pickle_file:
            pickle.dump(to_pickle, pickle_file)
    voltages = []
    for V in parameter_dict["voltage_sweep"]:
        voltages.append(V)
    proc_IDs = parameter_dict["proc_IDs"]
    jobs_list = [
        voltages[i : i + (int(np.ceil(len(voltages) / len(proc_IDs))))]
        for i in range(
            0, len(voltages), int(np.ceil(len(voltages) / float(len(proc_IDs))))
        )
    ]
    running_jobs = []
    output_dir = os.path.join(
        parameter_dict["output_device_dir"], parameter_dict["device_morphology"], "KMC"
    )
    print("Writing job pickles for each CPU...")
    for proc_ID, jobs in enumerate(jobs_list):
        pickle_name = os.path.join(output_dir, "KMC_data_{:02d}.pickle".format(proc_ID))
        with open(pickle_name, "wb+") as pickle_file:
            pickle.dump(jobs, pickle_file)
        print(
            "KMC jobs for proc_ID",
            proc_ID,
            "written to KMC_data_{:02d}.pickle".format(proc_ID),
        )
        # Open the required processes to execute the KMC jobs
        # Random seeding is a little weird here. If we don't generate a random
        # seed in the child process, it will just use the system time. So, we
        # generate a seed here to get the same random number stream each time,
        # and then feed the child process a new seed from the random number
        # stream. This way, we ensure that each child process has a different
        # random number stream to the other processes, but it's the same stream
        # every time we run the program.
        child_seed = np.random.randint(0, 2 ** 32)
        # Previous run command:
        run_command = [
            "python",
            SINGLE_RUN_DEVICE_KMC_FILE,
            output_dir,
            str(proc_ID),
            str(child_seed),
        ]
        print(run_command)
        running_jobs.append(sp.Popen(run_command))
    # Wait for all jobs to complete
    [p.wait() for p in running_jobs]
    print("All KMC jobs completed!")
    # Combine results if required.


if __name__ == "__main__":
    try:
        pickle_file = sys.argv[1]
    except:
        print(
            "Please specify the pickle file to load to continue the pipeline from"
            " this point."
        )
    _, _, _, parameter_dict, _ = hf.load_pickle(pickle_file)
    main(parameter_dict)
