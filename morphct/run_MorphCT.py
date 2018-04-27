import os
import shutil
import glob
import subprocess as sp

from morphct.code import helper_functions as hf
from morphct.code import fine_grainer
from morphct.code import run_HOOMD
from morphct.code import obtain_chromophores
from morphct.code import execute_ZINDO
from morphct.code import transfer_integrals
from morphct.code import mobility_KMC
from morphct.code import device_KMC


class simulation:
    def __init__(self, **kwargs):
        parameter_dict = {}
        # Read in all of the keyword arguments from the par file
        for key, value in kwargs.items():
            self.__dict__[key] = value
        # Obtain the slurm job ID (if there is one)
        self.slurm_job_ID = self.get_slurm_ID()
        # Parse the parameter file to get more useful file locations
        if self.morphology is not None:
            self.input_morphology_file = self.input_morph_dir + '/' + self.morphology
            self.output_morphology_directory = self.output_morph_dir + '/' + self.morphology[:-4]
        if self.device_morphology is not None:
            self.input_device_file = self.input_device_dir + '/' + self.device_morphology
            self.output_device_directory = self.output_device_dir + '/' + self.device_morphology
        # Add all the parameters to the parameterDict, which will be used to
        # send everything between classes
        for key, value in self.__dict__.items():
            if key in ['os', 'sys']:
                continue
            parameter_dict[key] = value
        # Make the correct directory tree
        self.make_dir_tree()
        if self.morphology is not None:
            # Copy the current code and the parameter file for safekeeping
            self.copy_code()
            if self.execute_fine_graining is False:
                # Load any previous data to allow us to run individual phases
                try:
                    pickle_data = hf.load_pickle(self.output_morphology_directory + '/code/'
                                                 + self.morphology[:-4] + '.pickle')
                    AA_morphology_dict = pickle_data[0]
                    CG_morphology_dict = pickle_data[1]
                    CG_to_AAID_master = pickle_data[2]
                    previous_parameter_dict = pickle_data[3]
                    chromophore_list = pickle_data[4]
                    # Load in any parameters from the previous_parameter_dict
                    # that have not been already defined in the new
                    # parameter_dict (e.g. CG_type_mappings):
                    for key, previous_value in previous_parameter_dict.items():
                        if key not in list(parameter_dict.keys()):
                            parameter_dict[key] = previous_value
                except:
                    print("PICKLE NOT FOUND, EXECUTING FINE GRAINING TO OBTAIN REQUIRED PARAMETERS...")
                    self.execute_fine_graining = False
                    returned_data = fine_grainer.morphology(self.input_morphology_file, self.morphology[:-4],
                                                            parameter_dict, []).analyse_morphology()
                    AA_morphology_dict = returned_data[0]
                    CG_morphology_dict = returned_data[1]
                    CG_to_AAID_master = returned_data[2]
                    parameter_dict = returned_data[3]
                    chromophore_list = returned_data[4]
            # Now begin running the code based on user's flags
            else:
                print("---=== BACKMAPPING COARSE-GRAINED SITES... ===---")
                returned_data = fine_grainer.morphology(self.input_morphology_file, self.morphology[:-4],
                                                        parameter_dict, []).analyse_morphology()
                AA_morphology_dict = returned_data[0]
                CG_morphology_dict = returned_data[1]
                CG_to_AAID_master = returned_data[2]
                parameter_dict = returned_data[3]
                chromophore_list = returned_data[4]
                print("---=== BACKMAPPING COMPLETED ===---")
            if self.execute_molecular_dynamics is True:
                print("---=== EQUILIBRATING FINE-GRAINED MORPHOLOGY... ===---")
                returned_data = run_HOOMD.main(AA_morphology_dict, CG_morphology_dict, CG_to_AAID_master,
                                               parameter_dict, chromophore_list)
                AA_morphology_dict = returned_data[0]
                CG_morphology_dict = returned_data[1]
                CG_to_AAID_master = returned_data[2]
                parameter_dict = returned_data[3]
                chromophore_list = returned_data[4]
                print("---=== EQUILIBRATION COMPLETED ===---")
            if self.execute_obtain_chromophores is True:
                print("---=== IDENTIFYING CHROMOPHORES OF CHARGE CARRIER DELOCALISATION... ===---")
                returned_data = obtain_chromophores.main(AA_morphology_dict, CG_morphology_dict, CG_to_AAID_master,
                                                         parameter_dict, chromophore_list)
                AA_morphology_dict = returned_data[0]
                CG_morphology_dict = returned_data[1]
                CG_to_AAID_master = returned_data[2]
                parameter_dict = returned_data[3]
                chromophore_list = returned_data[4]
                print("---=== IDENTIFICATION COMPLETED ===---")
            if self.execute_zindo is True:
                print("---=== PERFORMING SEMI-EMPIRICAL ZINDO/S CALCULATIONS... ===---")
                returned_data = execute_ZINDO.main(AA_morphology_dict, CG_morphology_dict, CG_to_AAID_master,
                                                   parameter_dict, chromophore_list)
                AA_morphology_dict = returned_data[0]
                CG_morphology_dict = returned_data[1]
                CG_to_AAID_master = returned_data[2]
                parameter_dict = returned_data[3]
                chromophore_list = returned_data[4]
                print("---=== CALCULATIONS COMPLETED ===---")
            if self.execute_calculate_transfer_integrals is True:
                print("---=== DETERMINING ELECTRONIC TRANSFER INTEGRALS... ===---")
                returned_data = transfer_integrals.main(AA_morphology_dict, CG_morphology_dict, CG_to_AAID_master,
                                                        parameter_dict, chromophore_list)
                AA_morphology_dict = returned_data[0]
                CG_morphology_dict = returned_data[1]
                CG_to_AAID_master = returned_data[2]
                parameter_dict = returned_data[3]
                chromophore_list = returned_data[4]
                print("---=== DETERMINATION COMPLETED ===---")
            if self.execute_calculate_mobility is True:
                print("---=== EXECUTING KINETIC MONTE CARLO MOBILITY SIMULATIONS... ===---")
                returned_data = mobility_KMC.main(AA_morphology_dict, CG_morphology_dict, CG_to_AAID_master,
                                                  parameter_dict, chromophore_list)
                AA_morphology_dict = returned_data[0]
                CG_morphology_dict = returned_data[1]
                CG_to_AAID_master = returned_data[2]
                parameter_dict = returned_data[3]
                chromophore_list = returned_data[4]
                print("---=== EXECUTION COMPLETED ===---")
        else:
            # NEED TO PUT A CHECK IN HERE TO ENSURE THAT WE LOAD THE CORRECT MOLECULAR DATA IN
            if self.main_device_simulation is True:
                print("---=== EXECUTING KINETIC MONTE CARLO DEVICE SIMULATIONS... ===---")
                device_KMC.main(parameter_dict)
                print("---=== EXECUTION COMPLETED ===---")
        exit()

    def get_slurm_ID(self):
        # Use Squeue to determine the current slurm job number
        try:
            squeue_command = sp.Popen(['squeue', '-u', str(os.getenv('user')), '--sort=t'],
                                      stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
        # If Slurm is not installed...
        except OSError:
            return None
        # ...or if the squeue errors out, then return no slurm job ID
        if len(squeue_command[1]) != 0:
            print("StdErr not empty:", squeue_command[1])
            return None
        output_lines = squeue_command[0].decode().split('\n')
        # If the command ran, the output is sorted by ascending runtime, so
        # this job will be the most recent submission from the current user
        # which is outputLines[1]
        for element in output_lines[1].split(' '):
            if len(element) != 0:
                # First element come across is the jobID
                return int(element)

    def make_dir_tree(self):
        print("Sorting out directory structure...")
        # Delete any previous data if the user asked to
        # if self.overwriteCurrentData == True:
        #     sp.Popen('echo rm -rf '+self.output_directory+'/*', shell=True)
        #     # Make sure that the rm command has finished before moving on
        #     sp.Popen('rm -rf '+self.output_directory+'/*', shell=True).communicate()
        # Then, make sure that all the required directories are in place
        if self.morphology is not None:
            for directory_to_make in ['chromophores/input_orca/single', 'chromophores/input_orca/pair',
                                      'chromophores/output_orca/single', 'chromophores/output_orca/pair',
                                      'KMC', 'molecules', 'morphology', 'code']:
                print('mkdir -p ' + self.output_morphology_directory + '/' + directory_to_make)
                # Make sure that the mkdir command has finished before moving on
                os.makedirs(self.output_morphology_directory + '/' + directory_to_make, exist_ok=True)
        elif self.device_morphology is not None:
            if self.overwrite_current_data is True:
                print('rm -r ' + self.outputDeviceDirectory + '/')
                shutil.rmtree(self.output_device_directory + '/', ignore_errors=True)
            for device_directory_to_make in ['code', 'KMC', 'figures']:
                if device_directory_to_make == 'figures':
                    for potential_val in self.voltage_sweep:
                        directory = device_directory_to_make + '/' + str(potential_val)
                        print('mkdir -p ' + self.output_device_directory + '/' + directory)
                        os.makedirs(self.output_device_directory + '/' + directory, exist_ok=True)
                else:
                    print('mkdir -p ' + self.output_device_directory + '/' + device_directory_to_make)
                    os.makedirs(self.output_device_directory + '/' + device_directory_to_make, exist_ok=True)

    def copy_code(self):
        print("Copying code...")
        code_dir = os.getcwd() + '/code'
        if self.morphology is not None:
            print('cp ' + code_dir + '/*.py ' + self.output_morphology_directory + '/code/')
            print('cp ' + os.getcwd() + '/' + self.parameter_file + ' '
                  + self.output_morphology_directory + '/code/')
            print('cp ' + self.input_morphology_file + ' ' + self.output_morphology_directory + '/code/input.xml')
            shutil.copy(os.getcwd() + '/' + self.parameter_file, self.output_morphology_directory + '/code')
            for file_name in glob.glob(code_dir + '/*.py'):
                shutil.copy(file_name, self.output_morphology_directory + '/code/')
            shutil.copy(os.getcwd() + '/' + self.parameter_file, self.output_morphology_directory + '/code/')
            shutil.copy(self.input_morphology_file, self.output_morphology_directory + '/code/input.xml')
        elif self.device_morphology is not None:
            print('cp ' + code_dir + '/*.py ' + self.output_device_directory + '/code/')
            print('cp ' + os.getcwd() + '/' + self.parameter_file + ' ' + self.output_device_directory + '/code/')
            shutil.copy(os.getcwd() + '/' + self.parameter_file, self.output_device_directory + '/code')
            for file_name in glob.glob(code_dir + '/*.py'):
                shutil.copy(file_name, self.output_device_directory + '/code/')
            shutil.copy(os.getcwd() + '/' + self.parameter_file, self.output_device_directory + '/code/')
