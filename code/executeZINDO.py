import os
import sys
import numpy as np
import time as T
import helperFunctions
import subprocess as sp



def createInputFiles(chromophoreList, AAMorphologyDict, parameterDict):
    # Singles first
    for chromophore in chromophoreList:
        # Include the molecule terminating units on the required atoms of the chromophore
        terminatingGroupPositions = terminateMonomers(chromophore, parameterDict, AAMorphologyDict)
        writeOrcaInp(AAMorphologyDict, chromophore.AAIDs, [[0, 0, 0]] * len(chromophore.AAIDs), terminatingGroupPositions, [[0, 0, 0]] * len(terminatingGroupPositions), chromophore.orcaInput)
    # Then pairs of all neighbours
    for chromophore1 in chromophoreList:
        neighboursID = [neighbour[0] for neighbour in chromophore1.neighbours]
        neighboursImage = [neighbour[1] for neighbour in chromophore1.neighbours]
        for chromophore2 in chromophoreList:
            if chromophore2.ID not in neighboursID:
                continue
            # Find the correct relative image for the neighbour chromophore
            chromophore2Image = neighboursImage[neighboursID.index(chromophore2.ID)]
            # Find the dimer AAIDs and relative images for each atom
            AAIDs = chromophore1.AAIDs + chromophore2.AAIDs
            images = [[0, 0, 0] for i in range(len(chromophore1.AAIDs))] + [chromophore2Image for i in range(len(chromophore2.AAIDs))]
            # Now add the terminating groups to both chromophores
            terminatingGroupPositions1 = terminateMonomers(chromophore1, parameterDict, AAMorphologyDict)
            terminatingGroupImages1 = [[0, 0, 0] for i in range(len(terminatingGroupPositions1))]
            terminatingGroupPositions2 = terminateMonomers(chromophore2, parameterDict, AAMorphologyDict)
            terminatingGroupImages2 = [chromophore2Image for i in range(len(terminatingGroupPositions2))]
            # Update the ORCA input name
            inputName = chromophore1.orcaInput.replace('.inp', '-%04d.inp' % (chromophore2.ID))
            # Write the dimer input file
            writeOrcaInp(AAMorphologyDict, AAIDs, images, terminatingGroupPositions1 + terminatingGroupPositions2, terminatingGroupImages1 + terminatingGroupImages2, inputName)


def writeOrcaInp(AAMorphologyDict, AAIDs, images, terminatingGroupPosns, terminatingGroupImages, inputName):
    linesToWrite = []
    allAtomTypes = []
    allPositions = []
    # Format the atom positions ready for ORCA
    for index, atomID in enumerate(AAIDs):
        # Cut the integer bit off the atomType
        allAtomTypes.append(''.join([i for i in AAMorphologyDict['type'][atomID] if not i.isdigit()]))
        # Add in the correct periodic images to the position
        allPositions.append(AAMorphologyDict['unwrapped_position'][atomID] + np.array([(images[index][i] * [AAMorphologyDict['lx'], AAMorphologyDict['ly'], AAMorphologyDict['lz']][i]) for i in range(3)]))
    # Now add in the terminating Hydrogens
    for index, position in enumerate(terminatingGroupPosns):
        # Cut the integer bit off the atomType
        allAtomTypes.append('H')
        # Add in the correct periodic images to the position
        allPositions.append(position + np.array([(terminatingGroupImages[index][i] * [AAMorphologyDict['lx'], AAMorphologyDict['ly'], AAMorphologyDict['lz']][i]) for i in range(3)]))
    # Now geometrically centralize all of the atoms that are to be included in this input file to make it easier on ORCA
    centralPosition = np.array([np.average(np.array(allPositions)[:,0]), np.average(np.array(allPositions)[:,1]), np.average(np.array(allPositions)[:,2])])
    # Create the lines to be written in the input file
    for index, position in enumerate(allPositions):
        linesToWrite.append(" %s  %.5f  %.5f  %.5f\n" % (allAtomTypes[index], position[0] - centralPosition[0], position[1] - centralPosition[1], position[2] - centralPosition[2]))
    # Load the ORCA input template
    with open(os.getcwd() + '/templates/template.inp', 'r') as templateFile:
        inpFileLines = templateFile.readlines()
    # Insert the linesToWrite
    inpFileLines[-1:-1] = linesToWrite
    # Write the ORCA input file
    with open(inputName, 'w+') as orcaFile:
        orcaFile.writelines(inpFileLines)
    print "\rOrca Input File written as", inputName[helperFunctions.findIndex(inputName, '/')[-1]+1:],



def terminateMonomers(chromophore, parameterDict, AAMorphologyDict):
    # Get the connections to the terminating groups from the parXX.py
    terminatingBonds = [bond for bond in parameterDict['moleculeTerminatingConnections']]
    # Remove any termination connections that already exist (i.e. terminating unit at the end of the molecule)
    popList = []
    for bondNo, bond in enumerate(terminatingBonds):
        if bond[0] in np.array(np.array(chromophore.bonds)[:,0]):
            popList.append(bondNo)
    for index in sorted(popList, reverse = True):
        terminatingBonds.pop(index)
    # Because we didn't reorder the AAID list at any point, the integer in terminatingBond[1] should correspond
    # to the correct AAID in chromo.AAIDs
    AAIDsToAttachTo = [chromophore.AAIDs[index] for index in map(int, list(np.array(terminatingBonds)[:,1]))]

    #print chromophore.bonds
    #print AAIDsToAttachTo
    #print terminatingBonds
    #for i in AAIDsToAttachTo:
    #    print AAMorphologyDict['type'][i]

    # Now work out the positions of any bonded atoms for each of these terminating atoms to work out where we
    # should put the hydrogen
    newHydrogenPositions = []
    for terminatingAtomID in AAIDsToAttachTo:
        # To do this, find the relative positions of the bonded atoms...
        thisAtomPosition = AAMorphologyDict['unwrapped_position'][terminatingAtomID]
        averagePositionOfBondedAtoms = np.array([0.0, 0.0, 0.0])
        for bond in chromophore.bonds:
            if (bond[1] == terminatingAtomID):
                averagePositionOfBondedAtoms += np.array(AAMorphologyDict['unwrapped_position'][bond[2]]) - np.array(thisAtomPosition)
            elif (bond[2] == terminatingAtomID):
                averagePositionOfBondedAtoms += np.array(AAMorphologyDict['unwrapped_position'][bond[1]]) - np.array(thisAtomPosition)
        # ... then reverse that vector and whack on a hydrogen 1.06 A away in that direction
        newHydrogenPositions.append(thisAtomPosition + (- averagePositionOfBondedAtoms / np.linalg.norm(averagePositionOfBondedAtoms) * 1.06))
    # Return terminatingGroups (positions of those hydrogens to be added to the ORCA input)
    return newHydrogenPositions


def execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList):
    createInputFiles(chromophoreList, AAMorphologyDict, parameterDict)


if __name__ == "__main__":
    try:
        pickleFile = sys.argv[1]
    except:
        print "Please specify the pickle file to load to continue the pipeline from this point."
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList = helperFunctions.loadPickle(pickleFile)
    execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList)

















#
#
#def countOutputFiles(directory):
#    singleOutputs = os.listdir(directory+'/single/')
#    pairOutputs = os.listdir(directory+'/pair/')
#    orcaOutputs = 0
#    for fileName in singleOutputs:
#        if fileName[-4:] == '.out':
#            orcaOutputs += 1
#    for fileName in pairOutputs:
#        if fileName[-4:] == '.out':
#            orcaOutputs += 1
#    return orcaOutputs
#
#
#def execute(morphologyFile, slurmJobNumber):
#    morphologyName = morphologyFile[helperFunctions.findIndex(morphologyFile,'/')[-1]+1:]
#    inputDir = os.getcwd()+'/outputFiles/'+morphologyName+'/chromophores/inputORCA'
#    # Clear input files
#    try:
#        os.unlink(inputDir.replace('/inputORCA', '/*.log'))
#    except OSError:
#        pass
#    procIDs, jobsList = helperFunctions.getORCAJobs(inputDir)
#    numberOfInputs = sum([len(ORCAFilesToRun) for ORCAFilesToRun in jobsList])
#    print "Found", numberOfInputs, "ORCA files to run."
#    if numberOfInputs > 0:
#        # Create pickle file containing the jobs sorted by ProcID
#        pickleName = inputDir.replace('inputORCA', 'ORCAJobs.pickle')
#        with open(pickleName, 'w+') as pickleFile:
#            pickle.dump(jobsList, pickleFile)
#        print "ORCA job pickle written to", pickleName
#        if len(jobsList) <= len(procIDs):
#            procIDs = procIDs[:len(jobsList)]
#        runningJobs = []
#        for CPURank in procIDs:
#            print 'python '+os.getcwd()+'/code/singleCoreRunORCA.py '+os.getcwd()+'/outputFiles/'+morphologyName+' '+str(CPURank)+' &'
#            # os.system('python '+os.getcwd()+'/code/singleCoreRunORCA.py '+os.getcwd()+'/outputFiles/'+morphologyName+' '+str(CPURank)+' &')
#            runningJobs.append(sp.Popen(['python', str(os.getcwd())+'/code/singleCoreRunORCA.py', str(os.getcwd())+'/outputFiles/'+morphologyName, str(CPURank)]))
#        # Wait for all jobs to complete
#        exitCodes = [p.wait() for p in runningJobs]
#        os.system('rm '+inputDir.replace('inputORCA', 'ORCAJobs.pickle'))
#
#    # print "Checking for completed output files..."
#    # previousNumberOfOutputs = -1
#    # slurmCancel = False
#    # while True:
#    #     if slurmCancel == True:
#    #         print "Terminating program..."
#    #         os.system('rm '+inputDir.replace('inputORCA', 'ORCAJobs.pickle'))
#    #         os.system('scancel '+str(slurmJobNumber))
#    #         exit()
#    #     numberOfOutputs = countOutputFiles(os.getcwd()+'/outputFiles/'+morphologyName+'/chromophores/outputORCA')
#    #     if numberOfOutputs == numberOfInputs:
#    #         print "All", numberOfInputs, "output files present. Waiting one more iteration for current jobs to complete..."
#    #         slurmCancel = True
#    #     if numberOfOutputs == previousNumberOfOutputs:
#    #         print "No additional output files found this iteration - there are still", numberOfOutputs, "output files present. Is everything still working?"
#    #     previousNumberOfOutputs = numberOfOutputs
#    #     # Sleep for 20 minutes
#    #     T.sleep(1200)
#
#if __name__ == '__main__':
#    morphologyFile = sys.argv[1]
#    slurmJobNumber = sys.argv[2]
#    execute(morphologyFile, slurmJobNumber)
