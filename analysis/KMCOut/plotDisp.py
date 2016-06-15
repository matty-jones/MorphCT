import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import csv
import scipy.optimize
import scipy.stats


elementaryCharge = 1.60217657E-19 # C
kB = 1.3806488E-23 # m^{2} kg s^{-2} K^{-1}
temperature = 290 # K

def loadCSVs(CSVDir):
    CSVList = []
    completeCSVData = {}
    targetTimes = {}
    dataDirs = [CSVDir]
    # dataDirs = [CSVDir+'/attempt2']
    # for fileName in os.listdir(CSVDir):
    #     if 'attempt' in fileName:
    #         dataDirs.append(CSVDir+'/'+fileName)
    for dataDir in dataDirs:
        for fileName in os.listdir(dataDir):
            if ".csv" in fileName:
                CSVList.append(dataDir+'/'+fileName)
    for CSVFileName in CSVList:
        slashLocs = findIndex(CSVFileName, '/')
        CSVName = str(CSVFileName[slashLocs[-1]+1:])
        completeCSVData[CSVName] = []
        underscoreLoc = findIndex(CSVFileName, '_')
        dotLoc = findIndex(CSVFileName, '.')
        targetTimes[CSVName] = float(CSVFileName[underscoreLoc[0]+1:dotLoc[-1]])
        with open(CSVFileName, 'r') as CSVFile:
            CSVData = csv.reader(CSVFile, delimiter=',')
            for row in CSVData:
                completeCSVData[CSVName].append(map(float, row))
    targetTimesList = []
    CSVDataList = []
    for key in completeCSVData.keys():
        targetTimesList.append(targetTimes[key])
        CSVDataList.append(completeCSVData[key])
    return CSVDataList, targetTimesList
        
    
def plotHist(CSVFile, targetTime, CSVDir):
    # CSVArray = np.array(CSVFile)
    # MSD = np.average(CSVArray[:,1]**2)
    # meanTime = np.average(CSVArray[:,3])
    disps = []
    squaredDisps = []
    times = []
    for carrierNo, carrierData in enumerate(CSVFile):
        # if (carrierData[3] > targetTime*2) or (carrierData[3] < targetTime/2.0):
        #     # When we don't squeeze the HOMO distribution, some hops take an extraordinarily long time
        #     # which means that we end up with crazy long hop times. This check just makes sure that
        #     # we only take into account ones that we care about
        #     # NOTE: There is a negligible number of these for M01TI0
        #    continue
        if carrierData[2] != 1: #Skip single hops
            disps.append(carrierData[1])
            squaredDisps.append(carrierData[1]**2)
            times.append(carrierData[3])
    if len(squaredDisps) == 0:
        print CSVFile
        print CSVDir
        print targetTime
        raise SystemError('HALT')
    MSD = np.average(squaredDisps)
    meanTime = np.average(times)
    if len(disps) < 2:
        plt.clf()
        return None, None
    #plt.figure()
    plt.hist(disps, 20)
    plt.ylabel('Frequency')
    plt.xlabel('Displacement (m)')
    fileName = 'disp_%.2E.png' % (targetTime)
    plt.savefig(CSVDir+'/'+fileName)
    plt.clf()
    print "Figure saved as", CSVDir+"/"+fileName
    return MSD, meanTime



def linearFit(x, y):
    x = np.array(x)
    y = np.array(y)
    '''Fits a linear fit of the form y = mx + c to the data'''
    fitfunc = lambda params, x: params[0] * x    #create fitting function of form y = mx + c
    errfunc = lambda p, x, y: fitfunc(p, x) - y              #create error function for least squares fit

    init_a = 0.5                            #find initial value for a (gradient)
    init_p = np.array((init_a))  #bundle initial values in initial parameters
    #calculate best fitting parameters (i.e. m and b) using the error function
    p1, success = scipy.optimize.leastsq(errfunc, init_p.copy(), args = (x, y))
    f = fitfunc(p1, x)          #create a fit with those parameters
    print p1
    print f
    return p1, f


def plotMSD(times, MSDs, CSVDir):
    fit = np.polyfit(times, MSDs, 1)
    fitX = np.linspace(np.min(times), np.max(times), 100)
    gradient, intercept, rVal, pVal, stdErr = scipy.stats.linregress(times, MSDs)
    print "Fitting rVal =", rVal
    #gradient, F = linearFit(times, MSDs)
    #exit()
    fitY = (fitX*gradient)
    mobility = calcMobility(fitX, fitY)
    #plt.figure()
    plt.plot(times, MSDs)
    plt.plot(fitX, fitY, 'r')
    plt.xlabel('Time (s)')
    plt.ylabel('MSD (m'+r'$^{2}$)')
    plt.title('Mob = '+str(mobility)+' cm'+r'$^{2}$/Vs', y = 1.1)
    fileName = 'LinMSD.png'
    plt.savefig(CSVDir+'/'+fileName)
    plt.clf()
    print "Figure saved as", CSVDir+"/"+fileName
    plt.semilogx(times, MSDs)
    plt.semilogx(fitX, fitY, 'r')
    plt.xlabel('Time (s)')
    plt.ylabel('MSD (m'+r'$^{2}$)')
    plt.title('Mob = '+str(mobility)+' cm'+r'$^{2}$/Vs', y = 1.1)
    fileName = 'LogMSD.png'
    plt.savefig(CSVDir+'/'+fileName)
    plt.clf()
    print "Figure saved as", CSVDir+"/"+fileName
    return mobility


def calcMobility(linFitX, linFitY):
    diffusionCoeff = (linFitY[-1] - linFitY[0])/(linFitX[-1] - linFitX[0])
    # Use Einstein relation (include the factor of 1/6!! It is in the Carbone/Troisi 2014 paper)
    mobility = elementaryCharge*diffusionCoeff/(6*kB*temperature) # This is in m^{2} / Vs
    return mobility*(100**2)


def parallelSort(list1, list2):
    '''This function sorts a pair of lists by the first list in ascending order (for example, atom mass and corresponding position can be input, sorted by ascending mass, and the two lists output, where the mass[atom_i] still corresponds to position[atom_i]'''
    data = zip(list1, list2)
    data.sort()
    list1, list2 = map(lambda t: list(t), zip(*data))
    return list1, list2

    
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


if __name__ == "__main__":
    tempDirs = []
    for fileName in os.listdir(os.getcwd()):
        if fileName[0] == 'T':
            tempDirs.append(fileName)
    temps = []
    mobs = []
    plt.figure()
    for tempDir in tempDirs:
        CSVDir = os.getcwd()+'/'+tempDir
        completeCSVData, targetTimes = loadCSVs(CSVDir)
        times = []
        MSDs = []
        for index, CSVFile in enumerate(completeCSVData):
            MSD, meanTime = plotHist(CSVFile, targetTimes[index], CSVDir)
            if MSD == None:
                continue
            times.append(meanTime)
            MSDs.append(MSD)
        times, MSDs = parallelSort(times, MSDs)
        mobility = plotMSD(times, MSDs, CSVDir)
        print "---=== Mobility for this KMC run ===---"
        print "Mobility =", mobility, "cm^{2} / Vs"
        print "---=================================---"
        MLoc = findIndex(tempDir, 'M')
        temps.append(tempDir[1:MLoc[0]])
        mobs.append(mobility)

    plt.semilogy(temps, mobs)
    plt.xlabel('Temperature, Arb. U')
    plt.ylabel('Mobility, cm'+r'$^{2}$ '+'V'+r'$^{-1}$'+r's$^{-1}$')
    plt.savefig('./mobTemp.png')
    plt.close()
    print "Mobility curve saved as './mobTemp.png'"
