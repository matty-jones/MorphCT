import sys
import os
import matplotlib.pyplot as plt
import numpy as np

def plotEnergy(timesteps, yvals, mode, molNumber):
    plt.clf()
    if mode != 'Std':
        plt.plot(timesteps, yvals, 'ro')
        plt.ylabel(mode)
        plt.xlabel('Timestep')
        if 'Energy' in mode:
            fileName = mode[:-7]+'_'+molNumber+'.png'
        else:
            fileName = mode+'_'+molNumber+'.png'
    else:
        xvals = np.arange(len(yvals))
        plt.plot(xvals, yvals, 'bo')
        plt.ylabel('Standard Deviation in Total Energy[-10:]')
        plt.xlabel('Timestep')
        fileName = 'STD_'+molNumber+'.png'
    plt.savefig(fileName)

def loadLog(fileName):
    with open(fileName, 'r') as logFile:
        data = logFile.readlines()
    timestep = []
    PE = []
    KE = []
    pairE = []
    bondE = []
    angleE = []
    dihedralE = []
    temperature = []
    pressure = []
    volume = []
    for line in data[1:]:
        splitLine = line.split('\t')
        splitLine[-1] = splitLine[-1][:-1]
        timestep.append(int(splitLine[0]))
        PE.append(float(splitLine[1]))
        KE.append(float(splitLine[2]))
        pairE.append(float(splitLine[3]))
        bondE.append(float(splitLine[4]))
        angleE.append(float(splitLine[5]))
        dihedralE.append(float(splitLine[6]))
        try:
            temperature.append(float(splitLine[7]))
            pressure.append(float(splitLine[8]))
            volume.append(float(splitLine[9]))
        except IndexError:
            # Old log file
            temperature = [0]*len(timestep)
            pressure = [0]*len(timestep)
            volume = [0]*len(timestep)
    totalEnergy = []
    totalPotential = []
    totalEnergyStd = []
    for valueNo in range(len(PE)):
        totalEnergy.append(PE[valueNo] + KE[valueNo] + pairE[valueNo] + bondE[valueNo] + angleE[valueNo] + dihedralE[valueNo])
        totalPotential.append(PE[valueNo] + pairE[valueNo] + bondE[valueNo] + angleE[valueNo] + dihedralE[valueNo])
        if valueNo >= 9:
            totalEnergyStd.append(np.std(totalEnergy[valueNo-10:valueNo]))
    return timestep, PE, KE, pairE, bondE, angleE, dihedralE, temperature, pressure, volume, totalEnergy, totalPotential, totalEnergyStd

def findIndex(string, character):
    '''This function returns the locations of an inputted character in an inputted string'''
    index = 0
    locations = []
    while index < len(string):
        if string[index] == character:
            locations.append(index)
        index += 1
    if len(locations) == 0:
        return None
    return locations




        # # print "Maximum STD =", self.maxStandardDeviation, "Current STD =", currentStd, "Target STD =", 0.05*self.maxStandardDeviation
        # if currentStd > self.maxStandardDeviation:
        #     self.maxStandardDeviation = currentStd
        #     stdIncreasing = True
        #     self.consecutiveDumpPeriodsUnderTarget = 0
        # else:
        #     stdIncreasing = False
        # self.standardDeviation.append(currentStd)
        # if (stdIncreasing == False) and (len(self.standardDeviation) >= 10):
        #     if currentStd <= 0.05*self.maxStandardDeviation:
        #         self.consecutiveDumpPeriodsUnderTarget += 1
        #     else:
        #         self.consecutiveDumpPeriodsUnderTarget = 0
        #     if self.consecutiveDumpPeriodsUnderTarget == 10:
        #         raise ExitHoomd("Standard Deviation Condition Met", self.moleculeName)
        # return 0



if __name__ == "__main__":
    cwd = os.getcwd()
    logFileName = sys.argv[1]
    molName = logFileName[findIndex(logFileName,'_')[0]+1:-4]
    timestep, PE, KE, pairE, bondE, angleE, dihedralE, temperature, pressure, volume, totalE, totalPE, totalEStd = loadLog(cwd+'/'+logFileName)

    newHeuristic = []
    movingPointAverageDegree = 100
    for i in range(0, len(totalPE)-movingPointAverageDegree, movingPointAverageDegree):
        rollingSum = 0
        for j in range(0, movingPointAverageDegree, 1):
            rollingSum += totalPE[i+j]
        rollingSum /= movingPointAverageDegree
        newHeuristic.append(rollingSum)

    MPATimestep = list(np.arange(len(newHeuristic)))
    plotEnergy(MPATimestep[1:], newHeuristic[1:], 'MPA', molName)
    # exit()


    
    plotEnergy(timestep[10:-10], PE[10:-10], 'Potential Energy', molName)
    plotEnergy(timestep[10:-10], KE[10:-10], 'Kinetic Energy', molName)
    plotEnergy(timestep[10:-10], pairE[10:-10], 'Pair Energy', molName)
    plotEnergy(timestep[10:-10], bondE[10:-10], 'Bond Energy', molName)
    plotEnergy(timestep[10:-10], angleE[10:-10], 'Angle Energy', molName)
    plotEnergy(timestep[10:-10], dihedralE[10:-10], 'Dihedral Energy', molName)
    plotEnergy(timestep[10:-10], temperature[10:-10], 'Temperature', molName)
    plotEnergy(timestep[10:-10], pressure[10:-10], 'Pressure', molName)
    plotEnergy(timestep[10:-10], volume[10:-10], 'Volume', molName)
    plotEnergy(timestep[10:-10], totalE[10:-10], 'Total Energy', molName)
    plotEnergy(timestep[10:-10], totalPE[10:-10], 'Total_Potential Energy', molName)
    plotEnergy(timestep[10:-10], totalEStd[10:-10], 'Std', molName)
