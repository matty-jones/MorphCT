import os
import numpy as np
import multiprocessing as mp
import time as T
import helperFunctions

def execute(morphologyFile):
    morphologyName = morphologyFile[helperFunctions.findIndex(morphologyFile,'/')[-1]+1:]
    inputDir = os.getcwd()+'/outputFiles/'+morphologyName+'/chromophores/inputORCA'
    # Clear input files
    try:
        os.unlink(inputDir.replace('/inputORCA', '/*.log'))
    except OSError:
        pass
    procIDs, jobsList = helperFunctions.getORCAJobs(inputDir)
    print "Found", sum([len(ORCAFilesToRun) for ORCAFilesToRun in jobsList]), "ORCA files to run."
    print procIDs
    for CPURank in procIDs:
        print 'python '+os.getcwd()+'/code/singleCoreRunORCA.py '+os.getcwd()+'/outputFiles/'+morphologyName+' '+str(CPURank)+' &'
        os.system('python '+os.getcwd()+'/code/singleCoreRunORCA.py '+os.getcwd()+'/outputFiles/'+morphologyName+' '+str(CPURank)+' &')
    

if __name__ == '__main__':
    morphologyFile = sys.argv[1]
    execute(morphologyFile)
