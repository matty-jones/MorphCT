import numpy as np
import sys
import os
from morphct.code import helperFunctions
from morphct.definitions import SINGLE_ORCA_RUN_FILE
import csv
import copy
import subprocess as sp
import multiprocessing as mp
import pickle
import time as T
import glob


class ORCAError(Exception):
    def __init__(self, fileName):
        self.string = "No molecular orbital data present for "+str(fileName)
    def __str__(self):
        return self.string


def loadORCAOutput(fileName):
    with open(fileName, 'r') as orcaFile:
        dataFile = orcaFile.readlines()
    recordMOData = False
    orbitalData = []
    for line in dataFile:
        if 'ORBITAL ENERGIES' in line:
            # Next line begins the MO data
            recordMOData = True
            continue
        if recordMOData == True:
            if 'MOLECULAR ORBITALS' in line:
                # Don't need anything else from the output file
                break
            dataInLine = []
            for element in line.split(' '):
                if len(element) > 1:
                    try:
                        dataInLine.append(float(element))
                    except ValueError:
                        continue
            if len(dataInLine) == 4:
                orbitalData.append(dataInLine)
    for i in range(len(orbitalData)):
        if orbitalData[i][1] == 0:
            # This line is the first unoccupied orbital - i.e. LUMO
            LUMO = orbitalData[i][3]
            HOMO = orbitalData[i-1][3]
            HOMO_1 = orbitalData[i-2][3]
            LUMO_1 = orbitalData[i+1][3]
            # Don't need any other orbitals
            break
    if recordMOData == False:
        # Molecular orbital data not present in this file
        raise ORCAError(fileName)
    return HOMO_1, HOMO, LUMO, LUMO_1


def modifyORCAFiles(fileName, failedFile, failedCount, chromophoreList):
    if failedCount == 3:
        # Three lots of reruns without any successes, try to turn off SOSCF
        print(str(fileName)+": Three lots of reruns without any success - turning off SOSCF to see if that helps...")
        turnOffSOSCF(failedFile)
    elif failedCount == 6:
        # Still no joy - increase the number of SCF iterations and see if convergence was just slow
        print(str(fileName)+": Six lots of reruns without any success - increasing the number of SCF iterations to 500...")
        increaseIterations(failedFile)
    elif failedCount == 9:
        # Finally, turn down the SCF tolerance
        print(str(fileName)+": Nine lots of reruns without any success - decreasing SCF tolerance (sloppySCF)...")
        reduceTolerance(failedFile)
    elif failedCount == 12:
        print(str(fileName)+": Failed to rerun ORCA 12 times, one final thing that can be done is to change the numerical accuracy...")
        revertORCAFiles(failedFile)
        increaseGrid(failedFile)
    elif failedCount == 15:
        print(str(fileName)+": Failed to rerun ORCA 15 times. Will try high numerical accuracy with no SOSCF as a last-ditch effort...")
        increaseGridNoSOSCF(failedFile)
    elif failedCount == 18:
        # SERIOUS PROBLEM
        print(str(fileName)+": Failed to rerun ORCA 18 times, even with all the input file tweaks. Examine the geometry - it is most likely unreasonable.")
        fileString = fileName[[index for index, char in enumerate(fileName) if char == "/"][-1] + 1:-4]
        for chromoString in fileString.split('-'):
            chromoID = int(chromoString)
            print("AAIDs for chromophore", chromoID)
            print(chromophoreList[chromoID].AAIDs)
        print("Reverting "+str(fileName)+" back to its original state...")
        revertORCAFiles(failedFile)
        return 1
    return 0


def turnOffSOSCF(inputFile):
    with open(inputFile, 'r') as fileName:
        originalLines = fileName.readlines()
    originalLines[3] = '!ZINDO/S NoSOSCF\n'
    with open(inputFile, 'w+') as fileName:
        fileName.writelines(originalLines)


def reduceTolerance(inputFile):
    with open(inputFile, 'r') as fileName:
        originalLines = fileName.readlines()
    originalLines[3] = '!ZINDO/S NoSOSCF SloppySCF\n'
    with open(inputFile, 'w+') as fileName:
        fileName.writelines(originalLines)


def increaseIterations(inputFile):
    with open(inputFile, 'r') as fileName:
        originalLines = fileName.readlines()
    originalLines.append('\n%scf MaxIter 500 end')
    with open(inputFile, 'w+') as fileName:
        fileName.writelines(originalLines)


def increaseGrid(inputFile):
    with open(inputFile, 'r') as fileName:
        originalLines = fileName.readlines()
    originalLines[3] = '!ZINDO/S SlowConv Grid7 NoFinalGrid\n'
    with open(inputFile, 'w+') as fileName:
        fileName.writelines(originalLines)


def increaseGridNoSOSCF(inputFile):
    with open(inputFile, 'r') as fileName:
        originalLines = fileName.readlines()
    originalLines[3] = '!ZINDO/S SlowConv Grid7 NoFinalGrid NoSOSCF SloppySCF\n'
    originalLines.append('\n%scf MaxIter 500 end')
    with open(inputFile, 'w+') as fileName:
        fileName.writelines(originalLines)


def revertORCAFiles(inputFile):
    with open(inputFile, 'r') as fileName:
        originalLines = fileName.readlines()
    originalLines[3] = '! ZINDO/S\n'
    for lineNo in range(len(originalLines)):
        # REMOVE THE SCF ITER
        if "%scf MaxIter" in originalLines[lineNo]:
            originalLines.pop(lineNo)
            break
    with open(inputFile, 'w+') as fileName:
        fileName.writelines(originalLines)


def rerunFails(failedChromoFiles, parameterDict, chromophoreList):
    print("")
    print(failedChromoFiles)
    print("There were", len(list(failedChromoFiles.keys())), "failed jobs.")
    procIDs = parameterDict['procIDs']
    outputDir = parameterDict['outputMorphDir'] + '/' + parameterDict['morphology'][:-4]
    popList = []
    permanentlyFailed = {}
    # Firstly, modify the input files to see if numerical tweaks make ORCA happier
    for failedFile, failedData in failedChromoFiles.items():
        failedCount = failedData[0]
        errorCode = modifyORCAFiles(failedFile, outputDir + '/chromophores/inputORCA/' + failedFile.replace('.out', '.inp'), failedCount, chromophoreList)
        if errorCode == 1:
            # Don't delete the elements from the list here because we're still trying to iterate over this dict and it cannot change length!
            popList.append(failedFile)
            permanentlyFailed[failedFile] = failedData
    # Now pop the correct elements from the failedChromoFiles dict
    for failedFile in popList:
        failedChromoFiles.pop(failedFile)
    # If there are no files left, then everything has failed so this function has completed its task
    if len(failedChromoFiles) == 0:
        return failedChromoFiles
    # Otherwise, rerun those failed files.
    # First, find the correct locations of the input Files
    inputFiles = [outputDir + '/chromophores/inputORCA/' + fileName.replace('.out', '.inp') for fileName in list(failedChromoFiles.keys())]
    # As before, split the list of reruns based on the number of processors
    jobsList = [inputFiles[i:i + (int(np.ceil(len(inputFiles) / len(procIDs)))) + 1] for i in range(0, len(inputFiles), int(np.ceil(len(inputFiles)/float(len(procIDs)))))]
    print(jobsList)
    # Write the jobs pickle for singleCoreRunORCA to obtain
    with open(outputDir + '/chromophores/ORCAJobs.pickle', 'wb+') as pickleFile:
        pickle.dump(jobsList, pickleFile)
    # Now rerun ORCA
    if len(jobsList) <= len(procIDs):
        procIDs = procIDs[:len(jobsList)]
    runningJobs = []
    for CPURank in procIDs:
        runningJobs.append(sp.Popen(['python', SINGLE_ORCA_RUN_FILE, outputDir, str(CPURank), '1'])) # The final argument here tells ORCA to ignore the presence of the output file and recalculate
    # Wait for running jobs to finish
    [p.wait() for p in runningJobs]
    # Finally, return the failed files list to the main failure handler to see if we need to iterate
    return failedChromoFiles, permanentlyFailed


def calculateDeltaE(chromophoreList, chromo1ID, chromo2ID):
    chromo1 = chromophoreList[chromo1ID]
    chromo2 = chromophoreList[chromo2ID]
    #### NOTE: SANITY CHECK  ####
    if (chromo1.ID != chromo1ID) or (chromo2.ID != chromo2ID):
        print("chromo1.ID (" + str(chromo1.ID) + ") != chromo1ID (" + str(chromo1ID) + "), or chromo2.ID (" + str(chromo2.ID) + ") != chromo2ID (" + str(chromo2ID) + ")! CHECK CODE!")
        exit()
    #### END OF SANITY CHECK ####
    if chromo1.species == 'Donor':
        # Hole transporter
        chromo1E = chromo1.HOMO
    elif chromo1.species == 'Acceptor':
        # Electron transporter
        chromo1E = chromo1.LUMO
    if chromo2.species == 'Donor':
        # Hole transporter
        chromo2E = chromo2.HOMO
    elif chromo2.species == 'Acceptor':
        # Electron transporter
        chromo2E = chromo2.LUMO
    #### NOTE: SANITY CHECK  ####
    if chromo1.species != chromo2.species:
        print("chromo1.species (" + str(chromo1.species) + ") != chromo2.species (" + str(chromo2.species) + ")! CHECK CODE!")
        exit()
    #### END OF SANITY CHECK ####
    return chromo2E - chromo1E, chromo1.species


def calculateTI(orbitalSplitting, deltaE):
    # Use the energy splitting in dimer method to calculate the electronic transfer integral in eV
    if deltaE**2 > orbitalSplitting**2:
        # Avoid an imaginary TI by returning zero.
        # (Could use KOOPMAN'S APPROXIMATION here if desired)
        TI = 0
    else:
        TI = 0.5 * np.sqrt((orbitalSplitting**2) - (deltaE**2))
    return TI


def updateSingleChromophoreList(chromophoreList, parameterDict):
    orcaOutputDir = parameterDict['outputMorphDir'] + '/' + parameterDict['morphology'][:-4] + '/chromophores/outputORCA/'
    # NOTE: This can possibly be done by recursively iterating through the neighbourlist of each chromophore, but I
    # imagine Python will whinge about the levels of recursion, so for now I'll just go through every chromophore twice.
    # Firstly, set the energy levels for each single chromophore, rerunning them if they fail.
    failedSingleChromos = {}  # Has the form {'FileName': [failCount, locationInChromophoreList]}
    for chromoLocation, chromophore in enumerate(chromophoreList):
        fileName = 'single/%04d.out' % (chromophore.ID)
        print("\rDetermining energy levels for", fileName, end=' ')
        sys.stdout.flush()
        # Update the chromophores in the chromophoreList with their energyLevels
        try:
            chromophore.HOMO_1, chromophore.HOMO, chromophore.LUMO, chromophore.LUMO_1 = loadORCAOutput(orcaOutputDir + fileName)
        except ORCAError:
            failedSingleChromos[fileName] = [1, chromoLocation]
            continue
    print("")
    # Rerun any failed ORCA jobs
    while len(failedSingleChromos) > 0:
        failedSingleChromos, permanentlyFailed = rerunFails(failedSingleChromos, parameterDict, chromophoreList)
        if len(permanentlyFailed) > 0:
            print(permanentlyFailed)
            print("--== CRITICAL ERROR ==--")
            print("THE ABOVE SINGLE-CHROMOPHORE SYSTEMS FAILED PERMANENTLY. THESE NEED FIXING/REMOVING FROM THE SYSTEM BEFORE ANY FURTHER DATA CAN BE OBTAINED.")
            raise SystemError("TERMINATING...")
        successfulReruns = []
        # Now check all of the files to see if we can update the chromophoreList
        for chromoName, chromoData in failedSingleChromos.items():
            print("Checking previously failed", chromoName)
            chromoID = chromoData[1]
            try:
                # Update the chromophore data in the chromophoreList
                chromophoreList[chromoID].HOMO_1, chromophoreList[chromoID].HOMO, chromophoreList[chromoID].LUMO, chromophoreList[chromoID].LUMO_1 = loadORCAOutput(orcaOutputDir + chromoName)
                # This chromophore didn't fail, so remove it from the failed list
                successfulReruns.append(chromoName)
            except ORCAError:
                # This chromophore failed so increment its fail counter
                failedSingleChromos[chromoName][0] += 1
                continue
        for chromoName in successfulReruns:
            failedSingleChromos.pop(chromoName)
    print("")
    # Finally, delete any of the files that need to be deleted.
    if parameterDict['removeORCAInputs'] is True:
        print("Deleting ORCA input files...")
        for fileName in glob.glob(orcaOutputDir.replace('outputORCA', 'inputORCA') + 'single/*.*'):
            os.remove(fileName)
    if parameterDict['removeORCAOutputs'] is True:
        print("Deleting ORCA output files...")
        for fileName in glob.glob(orcaOutputDir + 'single/*.*'):
            os.remove(fileName)
    return chromophoreList


def updatePairChromophoreList(chromophoreList, parameterDict):
    # Now that all the single chromophore energy levels are done, iterate through again and check the neighbours,
    # rerunning the pair file if it failed (which it won't have done because all my chromophores are delicious now).
    orcaOutputDir = parameterDict['outputMorphDir'] + '/' + parameterDict['morphology'][:-4] + '/chromophores/outputORCA/'
    failedPairChromos = {}
    for chromoLocation, chromophore in enumerate(chromophoreList):
        neighbourIDs = [neighbourData[0] for neighbourData in chromophore.neighbours]
        for neighbourLoc, neighbourID in enumerate(neighbourIDs):
            if chromophore.ID > neighbourID:
                continue
            fileName = 'pair/%04d-%04d.out' % (chromophore.ID, neighbourID)
            print("\rDetermining energy levels for", fileName, end=' ')
            sys.stdout.flush()
            try:
                dimerHOMO_1, dimerHOMO, dimerLUMO, dimerLUMO_1 = loadORCAOutput(orcaOutputDir + fileName)
                # Calculate the deltaE between the two single chromophores
                deltaE, species = calculateDeltaE(chromophoreList, chromophore.ID, neighbourID)
                # Calculate the TI using the ESD method
                if species == 'Donor':
                    TI = calculateTI(dimerHOMO - dimerHOMO_1, deltaE)
                elif species == 'Acceptor':
                    TI = calculateTI(dimerLUMO - dimerLUMO_1, deltaE)
                # Get the location of the current chromophore.ID in the neighbour's neighbourList
                reverseLoc = [neighbourData[0] for neighbourData in chromophoreList[neighbourID].neighbours].index(chromophore.ID)
                # Update both the current chromophore and the neighbour (for the reverse hop)
                chromophore.neighboursDeltaE[neighbourLoc] = deltaE
                chromophoreList[neighbourID].neighboursDeltaE[reverseLoc] = - deltaE
                chromophore.neighboursTI[neighbourLoc] = TI
                chromophoreList[neighbourID].neighboursTI[reverseLoc] = TI
                # DEBUG ASSERTIONS
                # Check list index corresponds to chromophore ID
                assert(chromoLocation == chromophoreList[chromoLocation].ID)
                assert(chromoLocation == chromophore.ID)
                # Check the neighbourLoc and reverseLoc give the correct chromophoreIDs
                assert(chromophoreList[chromophore.ID].neighbours[neighbourLoc][0] == chromophoreList[neighbourID].ID)
                assert(chromophoreList[neighbourID].neighbours[reverseLoc][0] == chromophoreList[chromophore.ID].ID)
                # Check the chromophoreList has been updated after updating the chromophore instance
                assert(chromophoreList[chromophore.ID].neighboursTI[neighbourLoc] == chromophore.neighboursTI[neighbourLoc])
                # Check the TI of the forward and backward hops are the same
                assert(chromophoreList[chromophore.ID].neighboursTI[neighbourLoc] == chromophoreList[neighbourID].neighboursTI[reverseLoc])
                # Check the chromophoreList has been updated after updating the chromophore instance
                assert(chromophoreList[chromophore.ID].neighboursDeltaE[neighbourLoc] == chromophore.neighboursDeltaE[neighbourLoc])
                # Check the DeltaE of the forward and backward hops are *= -1
                assert(chromophoreList[chromophore.ID].neighboursDeltaE[neighbourLoc] == -chromophoreList[neighbourID].neighboursDeltaE[reverseLoc])
                # END DEBUG ASSERTIONS
            except ORCAError:
                failedPairChromos[fileName] = [1, chromoLocation, neighbourID]
    print("")
    while len(failedPairChromos) > 0:
        failedPairChromos, permanentlyFailed = rerunFails(failedPairChromos, parameterDict, chromophoreList)
        if len(permanentlyFailed) > 0:
            print("--== WARNING ==--")
            print("The above chromophore-pair systems failed permanently. Setting their transfer integrals to zero, preventing these hops from ever taking place in the KMC.")
            for fileName, chromoData in permanentlyFailed.items():
                chromo1ID = chromoData[1]
                chromo2ID = chromoData[2]
                TI = 0.0
                deltaE = 0.0
                # Get the location of the neighbour's ID in the current chromophores's neighbourList
                neighbourLoc = [neighbourData[0] for neighbourData in chromophoreList[chromo1ID].neighbours].index(chromo2ID)
                # Get the location of the current chromophore's ID in the neighbour's neighbourList
                reverseLoc = [neighbourData[0] for neighbourData in chromophoreList[chromo2ID].neighbours].index(chromo1ID)
                # Update both the current chromophore and the neighbour (for the reverse hop)
                chromophoreList[chromo1ID].neighboursDeltaE[neighbourLoc] = deltaE
                chromophoreList[chromo2ID].neighboursDeltaE[reverseLoc] = - deltaE
                chromophoreList[chromo1ID].neighboursTI[neighbourLoc] = TI
                chromophoreList[chromo2ID].neighboursTI[reverseLoc] = TI
        successfulReruns = []
        for fileName, chromoData in failedPairChromos.items():
            print("Checking previously failed", fileName)
            chromo1ID = chromoData[1]
            chromo2ID = chromoData[2]
            try:
                dimerHOMO_1, dimerHOMO, dimerLUMO, dimerLUMO_1 = loadORCAOutput(orcaOutputDir + fileName)
                # Calculate the deltaE between the two single chromophores
                deltaE, species = calculateDeltaE(chromophoreList, chromophore.ID, neighbourID)
                # Calculate the TI using the ESD method
                if species == 'Donor':
                    TI = calculateTI(dimerHOMO - dimerHOMO_1, deltaE)
                elif species == 'Acceptor':
                    TI = calculateTI(dimerLUMO - dimerLUMO_1, deltaE)
                # Get the location of the neighbour's ID in the current chromophores's neighbourList
                neighbourLoc = [neighbourData[0] for neighbourData in chromophoreList[chromo1ID].neighbours].index(chromo2ID)
                # Get the location of the current chromophore's ID in the neighbour's neighbourList
                reverseLoc = [neighbourData[0] for neighbourData in chromophoreList[chromo2ID].neighbours].index(chromo1ID)
                # Update both the current chromophore and the neighbour (for the reverse hop)
                chromophoreList[chromo1ID].neighboursDeltaE[neighbourLoc] = deltaE
                chromophoreList[chromo2ID].neighboursDeltaE[reverseLoc] = - deltaE
                chromophoreList[chromo1ID].neighboursTI[neighbourLoc] = TI
                chromophoreList[chromo2ID].neighboursTI[reverseLoc] = TI
                # This rerun was successful so remove this chromophore from the rerun list
                successfulReruns.append(fileName)
                print(fileName, "was successful!")
                # DEBUG ASSERTIONS
                # Check the neighbourLoc and reverseLoc give the correct chromophoreIDs
                assert(chromophoreList[chromo1ID].neighbours[neighbourLoc][0] == chromophoreList[chromo2ID].ID)
                assert(chromophoreList[chromo2ID].neighbours[reverseLoc][0] == chromophoreList[chromo1ID].ID)
                # Check the TI of the forward and backward hops are the same
                assert(chromophoreList[chromo1ID].neighboursTI[neighbourLoc] == chromophoreList[chromo2ID].neighboursTI[reverseLoc])
                # Check the DeltaE of the forward and backward hops are *= -1
                assert(chromophoreList[chromo1ID].neighboursDeltaE[neighbourLoc] == -chromophoreList[chromo2ID].neighboursDeltaE[reverseLoc])
                # END DEBUG ASSERTIONS
            except ORCAError:
                # This dimer failed so increment its fail counter
                failedPairChromos[fileName][0] += 1
                print(fileName, "still failed, incrementing counter")
                continue
        print(len(failedPairChromos))
        for fileName in successfulReruns:
            failedPairChromos.pop(fileName)
        print(len(failedPairChromos))
    print("")
    # Finally, delete any of the files that need to be deleted.
    if parameterDict['removeORCAInputs'] is True:
        print("Deleting ORCA input files...")
        for fileName in glob.glob(orcaOutputDir.replace('outputORCA', 'inputORCA') + 'pair/*.*'):
            os.remove(fileName)
    if parameterDict['removeORCAOutputs'] is True:
        print("Deleting ORCA output files...")
        for fileName in glob.glob(orcaOutputDir + 'pair/*.*'):
            os.remove(fileName)
    return chromophoreList


def scaleEnergies(chromophoreList, parameterDict):
    # Shorter chromophores have significantly deeper HOMOs because they are treated as small molecules instead of chain segments.
    # To rectify this, find the average energy level for each chromophore and then map that average to the literature value
    # First, get the energy level data
    donorLevels = []
    acceptorLevels = []
    for chromo in chromophoreList:
        if (chromo.species == 'Donor'):
            donorLevels.append(chromo.HOMO)
        elif (chromo.species == 'Acceptor'):
            acceptorLevels.append(chromo.LUMO)
    if len(donorLevels) > 0:
        avHOMO = np.average(donorLevels)
        stdHOMO = np.std(np.array(donorLevels))
        deltaEHOMO = 0.0
    if len(acceptorLevels) > 0:
        avLUMO = np.average(acceptorLevels)
        stdLUMO = np.std(np.array(acceptorLevels))
        deltaELUMO = 0.0
    # Then add the lateral shift to to the energy levels to put the mean in line with the literature value
    # This is justified because we treat each chromophore in exactly the same way. Any deviation between the average of the calculated MOs and the literature one is therefore a systematic error arising from the short chromophore lengths and the frequency of the terminating groups in order to perform the DFT calculations.
    # By shifting the mean back to the literature value, we are accounting for this systematic error.
    if (parameterDict['literatureHOMO'] is None) and (parameterDict['literatureLUMO'] is None):
        # No energy level scaling necessary, move on to the target DoS width
        pass
    else:
        if (parameterDict['literatureHOMO'] is not None) and (len(donorLevels) > 0):
            deltaEHOMO = parameterDict['literatureHOMO'] - avHOMO
            donorLevels = list(np.array(donorLevels) + np.array([deltaEHOMO] * len(donorLevels)))
            avHOMO = parameterDict['literatureHOMO']
        if (parameterDict['literatureLUMO'] is not None) and (len(acceptorLevels) > 0):
            deltaELUMO = parameterDict['literatureLUMO'] - avLUMO
            acceptorLevels = list(np.array(acceptorLevels) + np.array([deltaELUMO] * len(acceptorLevels)))
            avLUMO = parameterDict['literatureLUMO']
        for chromo in chromophoreList:
            if (chromo.species == 'Donor'):
                deltaE = deltaEHOMO
            elif (chromo.species == 'Acceptor'):
                deltaE = deltaELUMO
            chromo.HOMO_1 += deltaE
            chromo.HOMO += deltaE
            chromo.LUMO += deltaE
            chromo.LUMO_1 += deltaE
    # Now squeeze the DoS of the distribution to account for the noise in these ZINDO/S calculations
    # Check the current STD of the DoS for both the donor and the acceptor, and skip the calculation if the current
    # STD is smaller than the literature value
    # First check the donor DoS
    if (parameterDict['targetDoSSTDHOMO'] is None) or (len(donorLevels) == 0):
        squeezeHOMO = False
    elif (parameterDict['targetDoSSTDHOMO'] > stdHOMO):
        squeezeHOMO = False
    else:
        squeezeHOMO = True
    # Then check the acceptor DoS
    if (parameterDict['targetDoSSTDLUMO'] is None) or (len(acceptorLevels) == 0):
        squeezeLUMO = False
    elif (parameterDict['targetDoSSTDLUMO'] > stdLUMO):
        squeezeLUMO = False
    else:
        squeezeLUMO = True
    if squeezeHOMO is True:
        for chromo in chromophoreList:
            if chromo.species == 'Donor':
                # Determine how many sigmas away from the mean this datapoint is
                sigma = (chromo.HOMO - avHOMO) / float(stdHOMO)
                # Calculate the new deviation from the mean based on the target STD and sigma
                newDeviation = parameterDict['targetDoSSTDHOMO'] * sigma
                # Work out the change in energy to be applied to meet this target energy level
                deltaE = (avHOMO + newDeviation) - chromo.HOMO
            else:
                continue
            # Apply the energy level displacement
            chromo.HOMO_1 += deltaE
            chromo.HOMO += deltaE
            chromo.LUMO += deltaE
            chromo.LUMO_1 += deltaE
    if squeezeLUMO is True:
        for chromo in chromophoreList:
            if chromo.species == 'Acceptor':
                # Determine how many sigmas away from the mean this datapoint is
                sigma = (chromo.LUMO - avLUMO) / float(stdLUMO)
                # Calculate the new deviation from the mean based on the target STD and sigma
                newDeviation = parameterDict['targetDoSSTDLUMO'] * sigma
                # Work out the change in energy to be applied to meet this target energy level
                deltaE = (avLUMO + newDeviation) - chromo.LUMO
            else:
                continue
            # Apply the energy level displacement
            chromo.HOMO_1 += deltaE
            chromo.HOMO += deltaE
            chromo.LUMO += deltaE
            chromo.LUMO_1 += deltaE
    return chromophoreList


def execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList):
    pickleName = parameterDict['outputMorphDir'] + '/' + parameterDict['morphology'][:-4] + '/code/' + parameterDict['morphology'][:-4] + '.pickle'
    # First, check that we need to examine the single chromophores
    runSingles = False
    if parameterDict['overwriteCurrentData'] is False:
        # Only perform this check if the user hasn't already specified to overwrite the data (in which case it runs anyway)
        # Run all singles if any of the single's data is missing (i.e. the HOMO level should suffice because all energy levels are updated at the same time, so we don't need to check all of them individually)
        for chromophore in chromophoreList:
            if chromophore.HOMO is None:
                runSingles = True
    if (runSingles is True) or (parameterDict['overwriteCurrentData'] is True):
        print("Beginning analysis of single chromophores...")
        chromophoreList = updateSingleChromophoreList(chromophoreList, parameterDict)
        # Now include any scaling to narrow the DoS or modulate the mean to match the literature HOMO/LUMO levels (which helps to negate the effect of short chromophores with additional hydrogens/terminating groups
        print("Scaling energies...")
        chromophoreList = scaleEnergies(chromophoreList, parameterDict)
        print("Single chromophore calculations completed. Saving...")
        helperFunctions.writePickle((AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList), pickleName)
    else:
        print("All single chromophore calculations already performed. Skipping...")
    # Then, check the pairs
    runPairs = False
    if parameterDict['overwriteCurrentData'] is False:
        for chromophore in chromophoreList:
            # Just check the first neighbour for each chromophore
            for neighbour in chromophore.neighboursTI:
                if neighbour is None:
                    runPairs = True
                    break
    if (runPairs is True) or (parameterDict['overwriteCurrentData'] is True):
        print("Beginning analysis of chromophore pairs...")
        chromophoreList = updatePairChromophoreList(chromophoreList, parameterDict)
        # DEBUG Testing - you can remove these as the assertions in updatePairChromophoreList should
        # already cover them, however they are fast and will ensure that there are no errors in the
        # chromophoreList after calculating the Tij and DeltaEijs
        TijError = checkForwardBackwardHopTij(chromophoreList)
        DeltaEError = checkForwardBackwardHopEij(chromophoreList)
        if TijError or DeltaEError:
            raise SystemError("Assertions failed, please address in code.")
        # END OF DEBUG Testing
        print("Pair chromophore calculations completed. Saving...")
        helperFunctions.writePickle((AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList), pickleName)
    else:
        print("All pair chromophore calculations already performed. Skipping...")
    return AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList


def checkForwardBackwardHopTij(chromophoreList):
    # Check reverse lookup: Tij === Tji
    donorErrors = 0
    acceptorErrors = 0
    for chromo1 in chromophoreList:
        chromo1ID = chromo1.ID
        for neighbourIndex, chromo2Deets in enumerate(chromo1.neighbours):
            chromo2ID = chromo2Deets[0]
            chromo1To2TI = chromo1.neighboursTI[neighbourIndex]
            # Sanity check
            assert(chromo2ID == chromophoreList[chromo2ID].ID)
            chromo2 = chromophoreList[chromo2ID]
            neighbour2Index = 0
            for neighbour2Index, chromo1Deets in enumerate(chromo2.neighbours):
                if chromo1Deets[0] != chromo1ID:
                    continue
                chromo2To1TI = chromo2.neighboursTI[neighbour2Index]
                break
            assert(chromo1.species == chromo2.species)
            try:
                assert(chromo1To2TI == chromo2To1TI)
                # Put other assertions in here
            except AssertionError:
                print("\n<ERROR FOUND>")
                print("1 to 2", chromo1To2TI)
                print("2 to 1", chromo2To1TI)
                print("Chromo 1 ID =", chromo1ID, "Chromo 2 ID =", chromo2ID)
                print("Chromo1 Neighbours: Look for index =", neighbourIndex, "in", chromo1.neighbours)
                print("Chromo2 Neighbours: Look for index =", neighbour2Index, "in", chromo2.neighbours)
                print("Chromo1 TIs: Look for index =", neighbourIndex, "in", chromo1.neighboursTI)
                print("Chromo2 TIs: Look for index =", neighbour2Index, "in", chromo2.neighboursTI)
                if chromo1.species == 'Donor':
                    donorErrors += 1
                elif chromo1.species == 'Acceptor':
                    acceptorErrors += 1
    if (donorErrors > 0) or (acceptorErrors > 0):
        print("--== CRITICAL ERROR ==--")
        print("\nThere were", donorErrors, "cases where Tij != Tji in the donor chromophores.")
        print("\nThere were", acceptorErrors, "cases where Tij != Tji in the acceptor chromophores.")
        return 1
    return 0


def checkForwardBackwardHopEij(chromophoreList):
    # Check reverse lookup: Delta Eij === -Delta Eji
    donorErrors = 0
    acceptorErrors = 0
    for chromophore in chromophoreList:
        chromoID = chromophore.ID
        for neighbourLoc, neighbourDeets in enumerate(chromophore.neighbours):
            neighbourID = neighbourDeets[0]
            assert(chromoID == chromophoreList[chromoID].ID)
            assert(neighbourID == chromophoreList[neighbourID].ID)
            # Get the location of the current chromophore.ID in the neighbour's neighbourList
            reverseLoc = [neighbourData[0] for neighbourData in chromophoreList[neighbourID].neighbours].index(chromophore.ID)
            assert(neighbourID == chromophoreList[chromoID].neighbours[neighbourLoc][0])
            assert(chromoID == chromophoreList[neighbourID].neighbours[reverseLoc][0])
            ## Update both the current chromophore and the neighbour (for the reverse hop)
            try:
                assert(chromophoreList[chromoID].neighboursDeltaE[neighbourLoc] == -chromophoreList[neighbourID].neighboursDeltaE[reverseLoc])
            except AssertionError:
                print("\nHOP FROM", chromoID, "TO", neighbourID)
                print(neighbourID, "should be here", chromophore.neighbours[neighbourLoc])
                print(chromoID, "should be here", chromophoreList[neighbourID].neighbours[reverseLoc])
                print("--== Transfer Integrals ==--")
                print("FORWARD:", chromophoreList[chromoID].neighboursTI[neighbourLoc], "BACKWARD:", chromophoreList[neighbourID].neighboursTI[reverseLoc])
                print("--== Delta Eij ==--")
                print("FORWARD:", chromophoreList[chromoID].neighboursDeltaE[neighbourLoc], "BACKWARD:", chromophoreList[neighbourID].neighboursDeltaE[reverseLoc])
                if chromophore.species == 'Donor':
                    donorErrors += 1
                elif chromophore.species == 'Acceptor':
                    acceptorErrors += 1
    if (donorErrors > 0) or (acceptorErrors > 0):
        print("--== CRITICAL ERROR ==--")
        print("\nThere were", donorErrors, "cases where Eij != -Eji in the donor chromophores.")
        print("\nThere were", acceptorErrors, "cases where Eij != -Eji in the acceptor chromophores.")
        return 1
    return 0


if __name__ == "__main__":
    try:
        pickleFile = sys.argv[1]
    except NameError:
        print("Please specify the pickle file to load to continue the pipeline from this point.")
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = helperFunctions.loadPickle(pickleFile)
    execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList)
