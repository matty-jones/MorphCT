import os
import sys
import shutil
import glob
import time as T
import subprocess as sp

sys.path.append(os.getcwd()+'/code')
import fineGrainer
import helperFunctions
try:
    import runHoomd
    import extractMol
except:
    print("HOOMD-Blue not found on this machine! As a result, runHoomd and extractMol will fail.")
import obtainChromophores
import executeZINDO
import transferIntegrals
import mobilityKMC
import deviceKMC


class simulation:
    def __init__(self, **kwargs):
        parameterDict = {}
        # Read in all of the keyword arguments from the par file
        for key, value in kwargs.items():
            self.__dict__[key] = value
        # Obtain the slurm job ID (if there is one)
        self.slurmJobID = self.getSlurmID()
        # Parse the parameter file to get more useful file locations
        if self.morphology is not None:
            self.inputMorphologyFile = self.inputMorphDir + '/' + self.morphology
            self.outputMorphologyDirectory = self.outputMorphDir + '/' + self.morphology[:-4]
        if self.deviceMorphology is not None:
            self.inputDeviceFile = self.inputDeviceDir + '/' + self.deviceMorphology
            self.outputDeviceDirectory = self.outputDeviceDir + '/' + self.deviceMorphology
        # Add all the parameters to the parameterDict, which will be used to send everything between classes
        for key, value in self.__dict__.items():
            if key in ['os', 'sys']:
                continue
            parameterDict[key] = value
        # Make the correct directory tree
        self.makeDirTree()
        if self.morphology is not None:
            # Copy the current code and the parameter file for safekeeping
            self.copyCode()
            if self.executeFinegraining is False:
            # Load any previous data to allow us to run individual phases
                try:
                    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, previousParameterDict, chromophoreList = helperFunctions.loadPickle(self.outputMorphologyDirectory+'/code/'+self.morphology[:-4]+'.pickle')
                    # Load in any parameters from the previousParameterDict that have not been already defined in the new parameterDict (e.g. CGTypeMappings):
                    for key, previousValue in previousParameterDict.items():
                        if key not in list(parameterDict.keys()):
                            parameterDict[key] = previousValue
                except:
                    print("PICKLE NOT FOUND, EXECUTING FINEGRAINING TO OBTAIN REQUIRED PARAMETERS...")
                    self.executeFinegraining = False
                    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = fineGrainer.morphology(self.inputMorphologyFile, self.morphology[:-4], parameterDict, []).analyseMorphology()
            # Now begin running the code based on user's flags
            else:
                print("---=== BACKMAPPING COARSE-GRAINED SITES... ===---")
                AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = fineGrainer.morphology(self.inputMorphologyFile, self.morphology[:-4], parameterDict, []).analyseMorphology()
                print("---=== BACKMAPPING COMPLETED ===---")
            if self.executeMolecularDynamics is True:
                print("---=== EQUILIBRATING FINE-GRAINED MORPHOLOGY... ===---")
                AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = runHoomd.execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList)
                print("---=== EQUILIBRATION COMPLETED ===---")
            if self.executeExtractMolecules is True:
                print("---=== EXTRACTING SINGLE MOLECULES FROM SYSTEM... ===---")
                AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = extractMol.execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList)
                print("---=== EXTRACTION COMPLETED ===---")
            if self.executeObtainChromophores is True:
                print("---=== IDENTIFYING CHROMOPHORES OF CHARGE CARRIER DELOCALISATION... ===---")
                AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = obtainChromophores.execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList)
                print("---=== IDENTIFICATION COMPLETED ===---")
            if self.executeZINDO is True:
                print("---=== PERFORMING SEMI-EMPIRICAL ZINDO/S CALCULATIONS... ===---")
                AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = executeZINDO.execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList)
                print("---=== CALCULATIONS COMPLETED ===---")
            if self.executeCalculateTransferIntegrals is True:
                print("---=== DETERMINING ELECTRONIC TRANSFER INTEGRALS... ===---")
                AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = transferIntegrals.execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList)
                print("---=== DETERMINATION COMPLETED ===---")
            if self.executeCalculateMobility is True:
                print("---=== EXECUTING KINETIC MONTE CARLO MOBILITY SIMULATIONS... ===---")
                AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = mobilityKMC.execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList)
                print("---=== EXECUTION COMPLETED ===---")
        else:
            # NEED TO PUT A CHECK IN HERE TO ENSURE THAT WE LOAD THE CORRECT MOLECULAR DATA IN
            if self.executeDeviceSimulation is True:
                print("---=== EXECUTING KINETIC MONTE CARLO DEVICE SIMULATIONS... ===---")
                deviceKMC.execute(parameterDict)
                print("---=== EXECUTION COMPLETED ===---")
        exit()

    def getSlurmID(self):
        # Use Squeue to determine the current slurm job number
        try:
            squeueCommand = sp.Popen(['squeue', '-u', os.getenv('USER'), '--sort=t'], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
        # If Slurm is not installed...
        except OSError:
            return None
        # ...or if the squeue errors out, then return no slurm job ID
        if len(squeueCommand[1]) != 0:
            print("StdErr not empty:", squeueCommand[1])
            return None
        outputLines = squeueCommand[0].decode().split('\n')
        # If the command ran, the output is sorted by ascending runtime, so this job will be the most recent submission from the current user which is outputLines[1]
        for element in outputLines[1].split(' '):
            if len(element) != 0:
                # First element come across is the jobID
                return int(element)

    def makeDirTree(self):
        print("Sorting out directory structure...")
        # Delete any previous data if the user asked to
        #if self.overwriteCurrentData == True:
        #    sp.Popen('echo rm -rf '+self.outputDirectory+'/*', shell=True)
        #    # Make sure that the rm command has finished before moving on
        #    sp.Popen('rm -rf '+self.outputDirectory+'/*', shell=True).communicate()
        # Then, make sure that all the required directories are in place
        # TODO: Remove the helperFunctions that mess around with the directory structure, do it all here instead.
        if self.morphology is not None:
            for directoryToMake in ['chromophores/inputORCA/single', 'chromophores/inputORCA/pair','chromophores/outputORCA/single', 'chromophores/outputORCA/pair', 'KMC', 'molecules', 'morphology', 'code']:
                print('mkdir -p ' + self.outputMorphologyDirectory + '/' + directoryToMake)
                # Make sure that the mkdir command has finished before moving on
                os.makedirs(self.outputMorphologyDirectory + '/' + directoryToMake, exist_ok=True)
        elif self.deviceMorphology is not None:
            if self.overwriteCurrentData is True:
                print('rm -r ' + self.outputDeviceDirectory + '/')
                shutil.rmtree(self.outputDeviceDirectory + '/', ignore_errors=True)
            for deviceDirectoryToMake in ['code', 'KMC', 'figures']:
                if deviceDirectoryToMake == 'figures':
                    for potentialVal in self.voltageSweep:
                        directory = deviceDirectoryToMake + '/' + str(potentialVal)
                        print('mkdir -p ' + self.outputDeviceDirectory + '/' + directory)
                        os.makedirs(self.outputDeviceDirectory + '/' + directory, exist_ok=True)
                else:
                    print('mkdir -p ' + self.outputDeviceDirectory + '/' + deviceDirectoryToMake)
                    os.makedirs(self.outputDeviceDirectory + '/' + deviceDirectoryToMake, exist_ok=True)


    def copyCode(self):
        print("Copying code...")
        codeDir = os.getcwd()+'/code'
        if self.morphology is not None:
            print('cp ' + codeDir + '/*.py ' + self.outputMorphologyDirectory + '/code/')
            print('cp ' + os.getcwd() + '/' + self.parameterFile + ' ' + self.outputMorphologyDirectory + '/code/')
            print('cp ' + self.inputMorphologyFile + ' ' + self.outputMorphologyDirectory + '/code/input.xml')
            shutil.copy(os.getcwd() + '/' + self.parameterFile, self.outputMorphologyDirectory + '/code')
            for fileName in glob.glob(codeDir + '/*.py'):
                shutil.copy(fileName, self.outputMorphologyDirectory+'/code/')
            shutil.copy(os.getcwd() + '/' + self.parameterFile, self.outputMorphologyDirectory + '/code/')
            shutil.copy(self.inputMorphologyFile, self.outputMorphologyDirectory + '/code/input.xml')
        elif self.deviceMorphology is not None:
            print('cp ' + codeDir + '/*.py ' + self.outputDeviceDirectory + '/code/')
            print('cp ' + os.getcwd() + '/' + self.parameterFile + ' ' + self.outputDeviceDirectory + '/code/')
            shutil.copy(os.getcwd() + '/' + self.parameterFile, self.outputDeviceDirectory + '/code')
            for fileName in glob.glob(codeDir + '/*.py'):
                shutil.copy(fileName, self.outputDeviceDirectory+'/code/')
            shutil.copy(os.getcwd() + '/' + self.parameterFile, self.outputDeviceDirectory + '/code/')
