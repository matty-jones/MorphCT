import sys
import os
import helperFunctions
import time as T

def writeToLogFile(logFile, string):
    with open(logFile, 'a+') as logWrite:
        logWrite.writelines(string+'\n')

if __name__ == '__main__':
    morphologyFile = sys.argv[1]
    CPURank = int(sys.argv[2])
    morphologyName = morphologyFile[helperFunctions.findIndex(morphologyFile, '/')[-1]+1:]
    inputDir = os.getcwd()+'/outputFiles/'+morphologyName+'/chromophores/inputORCA'
    logFile = inputDir.replace('/inputORCA', '/log_'+str(CPURank)+'.log')
    outputDir = os.getcwd()+'/outputFiles/'+morphologyName+'/chromophores/outputORCA'
    procIDs, jobsList = helperFunctions.getORCAJobs(inputDir)
    jobsToRun = jobsList[CPURank]
    writeToLogFile(logFile, 'Found '+str(len(jobsToRun))+' jobs to run.')
    t0 = T.time()
    for job in jobsToRun:
        t1 = T.time()
        writeToLogFile(logFile, 'Running job '+str(job)+'...')
        outputFileName = job.replace('.inp', '.out').replace('inputORCA', 'outputORCA')
        # Check if file exists already
        try:
            with open(outputFileName, 'r') as testFile:
                pass
            writeToLogFile(logFile, outputFileName+' already exists! Skipping...')
            continue
        except IOError:
            pass
        os.system(str(os.getcwd())+'/ORCA/orca '+str(job)+' > '+str(outputFileName))
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
        writeToLogFile(logFile, 'Job '+str(job)+' completed in '+elapsedTime+' '+timeunits)
        break
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
    writeToLogFile(logFile, 'All jobs completed in '+elapsedTime+' '+timeunits)
    writeToLogFile(logFile, 'Exiting normally...')
    
    
    
