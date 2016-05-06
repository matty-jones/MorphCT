import sys
import os
import multiprocessing as mp
import numpy as np
import subprocess as sp
import helperFunctions


def execute(morphologyFile):
    minTime = -11
    maxTime = -7
    morphologyName = morphologyFile[helperFunctions.findIndex(morphologyFile,'/')[-1]+1:]
    try:
        procIDs = list(np.arange(int(os.environ.get('SLURM_NPROCS'))))
    except (AttributeError, TypeError):
        procIDs = list(np.arange(mp.cpu_count()))
    timesList = np.logspace(minTime, maxTime, len(procIDs))
    for time in timesList:
        print 'python '+os.getcwd()+'/code/hoppingKMC.py '+os.getcwd()+'/outputFiles/'+morphologyName+' '+str(time)+' &'
        runningJobs.append(sp.Popen(['python', str(os.getcwd())+'/code/hoppingKMC.py', str(os.getcwd())+'/outputFiles/'+morphologyName, str(time)]))
    exitCodes = [p.wait() for p in runningJobs]
    

if __name__ == "__main__":
    morphologyFile = sys.argv[1]
    execute(morphologyFile)
