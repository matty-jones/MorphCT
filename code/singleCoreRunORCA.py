import sys
import os
import helperFunctions
import time as T
import subprocess as sp

def writeToFile(logFile, stringList, mode='logFile'):
    if mode == 'outputFile':
        openAs = 'w+'
    else:
        openAs = 'a+'
    with open(logFile, openAs) as logWrite:
        for line in stringList:
            logWrite.writelines(line+'\n')
        

if __name__ == '__main__':
    morphologyFile = sys.argv[1]
    CPURank = int(sys.argv[2])
    morphologyName = morphologyFile[helperFunctions.findIndex(morphologyFile, '/')[-1]+1:]
    orcaDir = os.getenv('ORCA_BIN', str(os.getcwd())+'/ORCA')
    orcaPath = orcaDir+'/orca'
    inputDir = os.getcwd()+'/outputFiles/'+morphologyName+'/chromophores/inputORCA'
    logFile = inputDir.replace('/inputORCA', '/ORCAlog_'+str(CPURank)+'.log')
    outputDir = os.getcwd()+'/outputFiles/'+morphologyName+'/chromophores/outputORCA'
    procIDs, jobsList = helperFunctions.getORCAJobs(inputDir)
    jobsToRun = jobsList[CPURank]
    writeToFile(logFile, ['Found '+str(len(jobsToRun))+' jobs to run.'])
    t0 = T.time()
    for job in jobsToRun:
        t1 = T.time()
        writeToFile(logFile, ['Running job '+str(job)+'...'])
        outputFileName = job.replace('.inp', '.out').replace('inputORCA', 'outputORCA')
        # Check if file exists already
        try:
            with open(outputFileName, 'r') as testFile:
                pass
            writeToFile(logFile, [outputFileName+' already exists! Skipping...'])
            continue
        except IOError:
            pass
        orcaJob = sp.Popen([str(orcaPath), str(job)], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
        jobPID = orcaJob.pid
        try:
            affinityJob = sp.Popen(['taskset', '-pc', str(CPURank), str(jobPID)], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
            #writeToFile(logFile, affinityJob[0].split('\n')) #stdOut for affinity set
            #writeToFile(logFile, affinityJob[1].split('\n')) #stdErr for affinity set
        except OSError:
            writeToFile(logFile, ["Taskset command not found, skipping setting of processor affinities..."])
        orcaShellOutput = orcaJob.communicate()
        # Write the outputFile:
        writeToFile(outputFileName, orcaShellOutput[0].split('\n'), mode = 'outputFile')
        #writeToFile(logFile, orcaShellOutput[0].split('\n')) #stdOut
        writeToFile(logFile, orcaShellOutput[1].split('\n')) #stdErr
        #os.system(orcaDir+'/orca '+str(job)+' > '+str(outputFileName))
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
        writeToFile(logFile, ['Job '+str(job)+' completed in '+elapsedTime+' '+timeunits])
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
    writeToFile(logFile, ['All jobs completed in '+elapsedTime+' '+timeunits])
    writeToFile(logFile, ['Exiting normally...'])
    
    
    
