import sys
import os
import multiprocessing as mp
import numpy as np
import subprocess as sp
import time as T

sys.path.append(os.getcwd()+'/code')
import helperFunctions


def getFilesList(direc):
    fileList = os.listdir(direc)
    morphologyFiles = []
    for fileName in fileList:
        if (fileName[-4:] == '.xml'):
            morphologyFiles.append(str(fileName))
    return morphologyFiles


def execute(morphologyFile):
    minTime = -11
    maxTime = -7
    morphologyName = morphologyFile[helperFunctions.findIndex(morphologyFile,'/')[-1]+1:]
    try:
        procIDs = list(np.arange(int(os.environ.get('SLURM_NPROCS'))))
    except (AttributeError, TypeError):
        procIDs = list(np.arange(mp.cpu_count()))
    timesList = np.logspace(minTime, maxTime, len(procIDs))
    runningJobs = []
    for time in timesList:
        print 'python '+os.getcwd()+'/code/hoppingKMC.py '+os.getcwd()+'/outputFiles/'+morphologyName+' '+str(time)+' &'
        runningJobs.append(sp.Popen(['python', str(os.getcwd())+'/code/hoppingKMC.py', str(os.getcwd())+'/outputFiles/'+morphologyName, str(time)]))
    exitCodes = [p.wait() for p in runningJobs]
    

if __name__ == "__main__":
    #### Put these into an input parameter file at some point along with everything else so that they are not defined here
    inputDir = 'inputCGMorphs'
    outputDir = 'outputFiles'
    currentDirContents = os.listdir(os.getcwd())
    if inputDir not in currentDirContents:
        os.makedirs(os.getcwd()+'/'+inputDir)
    if outputDir not in currentDirContents:
        os.makedirs(os.getcwd()+'/'+outputDir)
    inputDir = os.getcwd()+'/'+inputDir
    outputDir = os.getcwd()+'/'+outputDir
    
    morphologyFiles = sorted(getFilesList(inputDir))
    exitFlag = 0
    while exitFlag == 0:
        while True:
            print "\n---=== VALID MORPHOLOGY FILES ===---"
            for elementNo in range(len(morphologyFiles)):
                print str(elementNo)+"):", morphologyFiles[elementNo]
            print str(elementNo+1)+"): Exit program\n"
            # print "Valid files =", zip(datFiles, lammpstrjFiles)
            # AUTO RUNNING
            #runThisFile = raw_input("Please pick a file to run (integer, default = 0): ")
            try:
                runThisFile = str(sys.argv[1])
            except IndexError:
                runThisFile = '0'
            print "Selection =", runThisFile
            if len(runThisFile) == 0:
                runThisFile = 0
            else:
                try:
                    runThisFile = int(runThisFile)
                except:
                    print "Please enter an integer between 0 and", len(morphologyFiles)
                    continue
            if (runThisFile < 0) or (runThisFile > len(morphologyFiles)):
                print "Please enter an integer between 0 and", len(morphologyFiles)
                continue
            elif runThisFile == len(morphologyFiles):
                print "Exiting Program..."
                exitFlag = 1
                break
            break
        if exitFlag == 0:
            morphologyFile = outputDir+'/'+morphologyFiles[runThisFile][:-4]
            t0 = T.time()
            execute(morphologyFile)
            t1 = T.time()
            elapsedTime = float(t1) - float(t0)
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
            print "----------====================----------"
            print "KMC calculations completed in %.1f %s." % (float(elapsedTime), str(timeunits))
            print "----------====================----------"
            # Close program
            exitFlag = 1
            break
    print "Exitting program normally..."
