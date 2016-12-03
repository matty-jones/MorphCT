import sys
import signal
import traceback
import os
import numpy as np
import helperFunctions
import time as T
from scipy.sparse import lil_matrix
from scipy.sparse import find as findNonZero
import subprocess as sp
import cPickle as pickle


class terminationSignal:
    killSent = False
    def __init__(self):
        signal.signal(signal.SIGINT, self.catchKill)
        signal.signal(signal.SIGTERM, self.catchKill)

    def catchKill(self, signum, frame):
        self.killSent = True

class Terminate(Exception):
    '''This class is raised to terminate a KMC simulation'''
    def __init__(self, string):
        self.string = string

    def __str__(self):
        return self.string

def savePickle(jobsToRun, savePickleName):
    with open(savePickleName, 'w+') as pickleFile:
        pickle.dump(jobsToRun, pickleFile)
    helperFunctions.writeToFile(logFile, ['Pickle file saved successfully as ' + savePickleName + '!'])


def calculateDisplacement(initialPosition, finalPosition, finalImage, simDims):
    displacement = [0.0, 0.0, 0.0]
    for axis in range(3):
        displacement[axis] = (finalPosition[axis] - initialPosition[axis]) + (finalImage[axis] * (simDims[axis][1] - simDims[axis][0]))
    return np.linalg.norm(np.array(displacement))


if __name__ == '__main__':
    # NOTE: Remove print statements when done debugging. Everything should be written to the log file instead.




    KMCDirectory = sys.argv[1]
    CPURank = int(sys.argv[2])
    overwrite = False
    try:
        overwrite = bool(sys.argv[3])
    except:
        pass
    pickleFileName = KMCDirectory + '/carrierJobs.pickle'
    with open(pickleFileName, 'r') as pickleFile:
        jobsList = pickle.load(pickleFile)
    savePickleName = KMCDirectory + '/KMCData_' + str(CPURank) + '.pickle'
    logFile = KMCDirectory + '/KMClog_' + str(CPURank) + '.log'
    # Reset the log file
    with open(logFile, 'w+') as logFileHandle:
        pass
    jobsToRun = jobsList[CPURank]
    helperFunctions.writeToFile(logFile, ['Found ' + str(len(jobsToRun)) + ' jobs to run.'])
    # Set the affinities for this current process to make sure it's maximising available CPU usage
    currentPID = os.getpid()
    try:
        affinityJob = sp.Popen(['taskset', '-pc', str(CPURank), str(currentPID)], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
        # helperFunctions.writeToFile(logFile, affinityJob[0].split('\n')) #stdOut for affinity set
        # helperFunctions.writeToFile(logFile, affinityJob[1].split('\n')) #stdErr for affinity set
    except OSError:
        helperFunctions.writeToFile(logFile, ["Taskset command not found, skipping setting of processor affinity..."])
    # Now load the main morphology pickle (used as a workaround to obtain the chromophoreList without having to save it in each carrier [very memory inefficient!])
    pickleDir = KMCDirectory.replace('/KMC', '/code')
    for fileName in os.listdir(pickleDir):
        if 'pickle' in fileName:
            mainMorphologyPickleName = pickleDir + '/' + fileName
    helperFunctions.writeToFile(logFile, ['Found main morphology pickle file at ' + mainMorphologyPickleName + '! Loading data...'])
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList = helperFunctions.loadPickle(mainMorphologyPickleName)
    helperFunctions.writeToFile(logFile, ['Main morphology pickle loaded!'])
    # TEST THIS
    # Attempt to catch a kill signal to ensure that we save the pickle before termination
    killer = terminationSignal()
    # Save the pickle file as a `checkpoint' either every 100 carriers or every 1%, whichever is greater
    # While running (and in case the killer doesn't work, or if slurm sends a SIGKILL instead of a SIGTERM), save a checkpoint every 100 carriers or 1% of the number of total carriers to run, whichever is larger.
    checkpointCarriers = []
    onePercentOfJobs = int(len(jobsToRun) / 100)
    if onePercentOfJobs > 100:
        checkpointCarriers = list(map(int, np.arange(onePercentOfJobs, len(jobsToRun), onePercentOfJobs)))
    else:
        checkpointCarriers = list(map(int, np.arange(100, len(jobsToRun), 100)))
    t0 = T.time()
    try:
        for jobNumber, carrier in enumerate(jobsToRun):
            t1 = T.time()
            helperFunctions.writeToFile(logFile, ['Running carrier %d for %.2E...' % (carrier.ID, carrier.lifetime)])
            terminateSimulation = False
            while terminateSimulation is False:
                terminateSimulation = bool(carrier.calculateHop(chromophoreList))
                if killer.killSent is True:
                    raise Terminate('Kill command sent, terminating KMC simulation...')
            # Now the carrier has finished hopping, let's calculate its vitals
            initialPosition = carrier.initialChromophore.posn
            finalPosition = carrier.currentChromophore.posn
            finalImage = carrier.image
            simDims = carrier.simDims
            carrier.displacement = calculateDisplacement(initialPosition, finalPosition, finalImage, simDims)
            t2 = T.time()
            elapsedTime = float(t2) - float(t1)
            if elapsedTime < 60:
                timeunits = 'seconds.'
            elif elapsedTime < 3600:
                elapsedTime /= 60.0
                timeunits = 'minutes.'
            elif elapsedTime < 86400:
                elapsedTime /= 3600.0
                timeunits = 'hours.'
            else:
                elapsedTime /= 86400.0
                timeunits = 'days.'
            elapsedTime = '%.1f' % (float(elapsedTime))
            helperFunctions.writeToFile(logFile, ['Carrier hopped ' + str(carrier.noHops) + ' times into image ' + str(carrier.image) + ' for a displacement of ' + str(carrier.displacement) + ' in ' + str(elapsedTime) + ' ' + str(timeunits)])
            # Save the pickle file as a `checkpoint' either every 100 carriers or every 1%, whichever is greater
            if jobNumber in checkpointCarriers:
                print "Completed", jobNumber, "jobs. Making checkpoint at %3d%%" % (np.round(jobNumber + 1 / float(len(jobsToRun)) * 100))
                helperFunctions.writeToFile(logFile, ['Completed ' + str(jobNumber) + ' jobs. Making checkpoint at %3d%%' % (np.round(jobNumber / float(len(jobsToRun)) * 100))])
                savePickle(jobsToRun, savePickleName)
    except Exception as errorMessage:
        print traceback.format_exc()
        print "Saving the pickle file cleanly before termination..."
        helperFunctions.writeToFile(logFile, [str(errorMessage)])
        helperFunctions.writeToFile(logFile, ['Saving the pickle file cleanly before termination...'])
        savePickle(jobsToRun, savePickleName)
        exit()
    t3 = T.time()
    elapsedTime = float(t3) - float(t0)
    if elapsedTime < 60:
        timeunits = 'seconds.'
    elif elapsedTime < 3600:
        elapsedTime /= 60.0
        timeunits = 'minutes.'
    elif elapsedTime < 86400:
        elapsedTime /= 3600.0
        timeunits = 'hours.'
    else:
        elapsedTime /= 86400.0
        timeunits = 'days.'
    elapsedTime = '%.1f' % (float(elapsedTime))
    helperFunctions.writeToFile(logFile, ['All jobs completed in ' + elapsedTime + ' ' + timeunits])
    helperFunctions.writeToFile(logFile, ['Exiting normally...'])
