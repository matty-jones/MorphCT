import sys
import os
import helperFunctions
import time as T
import subprocess as sp
import pickle


if __name__ == '__main__':
    print(sys.argv)
    morphologyFile = sys.argv[1]
    CPURank = int(sys.argv[2])
    overwrite = False
    try:
        overwrite = bool(sys.argv[3])
    except:
        pass
    morphologyName = morphologyFile[helperFunctions.findIndex(morphologyFile, '/')[-1] + 1:]
    orcaDir = os.getenv('ORCA_BIN', str(os.getcwd()) + '/ORCA')
    orcaPath = orcaDir + '/orca'
    inputDir = os.getcwd() + '/outputFiles/' + morphologyName + '/chromophores/inputORCA'
    logFile = inputDir.replace('/inputORCA', '/ORCAlog_' + str(CPURank) + '.log')
    outputDir = os.getcwd() + '/outputFiles/' + morphologyName + '/chromophores/outputORCA'
    pickleFileName = inputDir.replace('inputORCA', 'ORCAJobs.pickle')
    with open(pickleFileName, 'r') as pickleFile:
        jobsList = pickle.load(pickleFile)
    jobsToRun = jobsList[CPURank]
    helperFunctions.writeToFile(logFile, ['Found ' + str(len(jobsToRun)) + ' jobs to run.'])
    t0 = T.time()
    for job in jobsToRun:
        t1 = T.time()
        helperFunctions.writeToFile(logFile, ['Running job ' + str(job) + '...'])
        outputFileName = job.replace('.inp', '.out').replace('inputORCA', 'outputORCA')
        # Check if file exists already
        if overwrite is False:
            try:
                with open(outputFileName, 'r') as testFile:
                    pass
                helperFunctions.writeToFile(logFile, [outputFileName + ' already exists! Skipping...'])
                continue
            except IOError:
                pass
        orcaJob = sp.Popen([str(orcaPath), str(job)], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
        jobPID = orcaJob.pid
        try:
            affinityJob = sp.Popen(['taskset', '-pc', str(CPURank), str(jobPID)], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
            # helperFunctions.writeToFile(logFile, affinityJob[0].split('\n')) #stdOut for affinity set
            # helperFunctions.writeToFile(logFile, affinityJob[1].split('\n')) #stdErr for affinity set
        except OSError:
            helperFunctions.writeToFile(logFile, ["Taskset command not found, skipping setting of processor affinities..."])
        orcaShellOutput = orcaJob.communicate()
        # Write the outputFile:
        helperFunctions.writeToFile(outputFileName, orcaShellOutput[0].split('\n'), mode='outputFile')
        # helperFunctions.writeToFile(logFile, orcaShellOutput[0].split('\n'))  # stdOut
        helperFunctions.writeToFile(logFile, orcaShellOutput[1].split('\n'))  # stdErr
        # os.system(orcaDir+'/orca '+str(job)+' > '+str(outputFileName))
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
        helperFunctions.writeToFile(logFile, ['Job ' + str(job) + ' completed in ' + elapsedTime + ' ' + timeunits])
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
