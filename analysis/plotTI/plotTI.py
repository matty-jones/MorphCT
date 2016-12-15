import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import csv
sys.path.append('../../code')
import helperFunctions

from scipy.optimize import curve_fit

elementaryCharge = 1.60217657E-19 # C
kB = 1.3806488E-23 # m^{2} kg s^{-2} K^{-1}
hbar = 1.05457173E-34 # m^{2} kg s^{-1}


def calculateLambdaij(chromoLength):
    # The equation for the internal reorganisation energy was obtained from the data given in
    # Johansson, E and Larsson, S; 2004, Synthetic Metals 144: 183-191.
    # External reorganisation energy obtained from 
    # Liu, T and Cheung, D. L. and Troisi, A; 2011, Phys. Chem. Chem. Phys. 13: 21461-21470
    lambdaExternal = 0.11 # eV
    if chromoLength < 12:
        lambdaInternal = 0.20826 - (chromoLength*0.01196)
    else:
        lambdaInternal = 0.06474
    lambdaeV = lambdaExternal+lambdaInternal
    return lambdaeV


def gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


def gaussFit(data):
    n = len(data)
    mean = np.mean(data)
    std = np.std(data)
    print "\n"
    print "Delta Eij stats: mean =", mean, "std =", std
    hist, binEdges = np.histogram(data, bins=100)
    fitArgs, fitConv = curve_fit(gaussian, binEdges[:-1], hist, p0=[1, mean, std])
    return binEdges, fitArgs


def plotHist(saveDir, yvals, mode, xvals=None, gaussBins=None, fitArgs=None):

    if mode == 'HOMO':
        plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('HOMO Level (eV)')
        plt.xlim([-5.5, -3.5])
        fileName = 'HOMODoS.pdf'

    elif mode == 'Bandgap':
        plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Bandgap (eV)')
        plt.xlim([3, 9])
        fileName = 'Bandgap.pdf'

    elif mode == 'BandgapLength':
        plt.scatter(xvals, yvals)
        plt.ylabel('Bandgap (eV)')
        plt.xlabel('Chromo Length (monomers)')
        fileName = 'BandgapLength.pdf'

    elif mode == 'Splitting':
        plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('HOMO Splitting (ev)')
        plt.xlim([0, 3.5])
        fileName = 'HOMOSplit.pdf'

    elif mode == 'TI':
        plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Transfer Integral (eV)')
        plt.xlim([0.0, 1.2])
        fileName = 'TI.pdf'

    elif mode == 'TITrimmed':
        plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Non-Zero Transfer Integral (eV)')
        plt.xlim([0.0, 1.2])
        fileName = 'TITrimmed.pdf'

    elif mode == 'Length':
        plt.scatter(xvals, yvals)
        plt.xlabel('HOMO Level (eV)')
        plt.ylabel('Chromo Length (monomers)')
        fileName = 'HOMOLength.pdf'

    elif mode == 'lambda':
        plt.scatter(xvals, yvals)
        plt.xlabel('Chromo Length (monomers)')
        plt.ylabel('Reorganisation energy (eV)')
        fileName = 'LambdaIJ.pdf'

    elif mode == 'deltaEij':
        n, bins, patches = plt.hist(yvals, 20)
        gaussY = gaussian(gaussBins[:-1], *fitArgs)
        scaleFactor = max(n)/max(gaussY)
        plt.plot(gaussBins[:-1], gaussY*scaleFactor, 'ro:')
        plt.ylabel('Frequency')
        plt.xlabel('Delta Eij (eV)')
        plt.xlim([-0.5, 0.5])
        fileName = 'deltaEij.pdf'

    elif mode == 'averageHOMO':
        plt.scatter(xvals, yvals)
        plt.xlabel('Chromo Length (monomers)')
        plt.ylabel('Average HOMO level (eV)')
        fileName = 'averageHOMO.pdf'

    elif mode == 'intraChainHop':
        if len(yvals) > 0:
            plt.hist(yvals, bins = np.logspace(1, 18, 40))
        plt.ylabel('Frequency')
        plt.xlabel('Intra-Chain Hop rate (s' + r'^{-1}' + ')')
        plt.gca().set_xscale('log')
        plt.xlim([1, 1E18])
        fileName = 'intrakij.pdf'

    elif mode == 'interChainHop':
        if len(yvals) > 0:
            plt.hist(yvals, bins = np.logspace(1, 18, 40))
        plt.ylabel('Frequency')
        plt.xlabel('Inter-Chain Hop rate (s' + r'$^{-1}$' + ')')
        plt.gca().set_xscale('log')
        plt.xlim([1, 1E18])
        fileName = 'interkij.pdf'

    elif mode == 'intraChainTI':
        if len(yvals) > 0:
            plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Intra-Chain TI (eV)')
        plt.xlim([0, 1.2])
        fileName = 'intraTij.pdf'

    elif mode == 'interChainTI':
        if len(yvals) > 0:
            plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Inter-Chain TI (eV)')
        plt.xlim([0, 1.2])
        fileName = 'interTij.pdf'

    elif mode == 'intraChainTITrim':
        if len(yvals) > 0:
            plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Intra-Chain TI (eV)')
        plt.xlim([0, 1.2])
        fileName = 'intraTijTrim.pdf'

    elif mode == 'interChainTITrim':
        if len(yvals) > 0:
            plt.hist(yvals, 20)
        plt.ylabel('Frequency')
        plt.xlabel('Inter-Chain TI (eV)')
        plt.xlim([0, 1.2])
        fileName = 'interTijTrim.pdf'

    elif mode == 'wobbey':
        plt.scatter(xvals, yvals)
        plt.xlabel('Chromo Length (monomers)')
        plt.ylabel('A (eV)')
        fileName = 'averageHOMO.pdf'

    elif mode == 'hopMix':
        plt.hist([yvals, xvals], bins = np.logspace(1, 18, 40), stacked = True, color = ['r', 'b'], label = ['Intra-chain', 'Inter-chain'])
        plt.ylabel('Frequency')
        plt.xlabel('Hopping Rate (s' + r'$^{-1}$' + ')')
        plt.xlim([1,1E18])
        plt.ylim([0,8000])
        plt.legend(loc = 2, prop = {'size':18})
        plt.gca().set_xscale('log')
        fileName = 'hoppingRateMixed.pdf'


    plt.savefig(saveDir+'/'+fileName)
    plt.clf()
    print "Figure saved as", saveDir+"/"+fileName


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


def calculateHopRate(lambdaij, Tij, deltaEij, T):
    # Semiclassical Marcus Hopping Rate Equation
    kij = ((2 * np.pi) / hbar) * (Tij ** 2) * np.sqrt(1.0 / (4 * lambdaij * np.pi * kB * T)) * np.exp(-((deltaEij + lambdaij)**2) / (4 * lambdaij * kB * T))
    return kij


if __name__ == "__main__":
    tempDirs = []
    for fileName in os.listdir(os.getcwd()):
        if fileName[0] == 'T':
            tempDirs.append(fileName)
    plt.figure()
    for tempDir in tempDirs:
        for fileName in os.listdir(os.getcwd() + '/' + tempDir):
            if "pickle" in fileName:
                mainMorphologyPickleName = os.getcwd() + '/' + tempDir + '/' + fileName
        AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList = helperFunctions.loadPickle(mainMorphologyPickleName)
        HOMOLevels = []
        LUMOLevels = []
        bandgap = []
        transferIntegrals = []
        transferIntegralsTrimmed = []
        deltaEij = []
        for chromophore in chromophoreList:
            if chromophore.species == 'Donor':
                HOMOLevels.append(chromophore.HOMO)
            else:
                LUMOLevels.append(chromophore.LUMO)
            bandgap.append(chromophore.LUMO - chromophore.HOMO)
            transferIntegrals += chromophore.neighboursTI
            for neighbourTI in chromophore.neighboursTI:
                if neighbourTI != 0:
                    transferIntegralsTrimmed.append(neighbourTI)
            deltaEij += chromophore.neighboursDeltaE

        binEdges, fitArgs = gaussFit(deltaEij)

        plotHist(tempDir, HOMOLevels, 'HOMO')
        plotHist(tempDir, transferIntegrals, 'TI')
        plotHist(tempDir, transferIntegralsTrimmed, 'TITrimmed')
        plotHist(tempDir, deltaEij, 'deltaEij', gaussBins=binEdges, fitArgs=fitArgs)
        plotHist(tempDir, bandgap, 'Bandgap')

        CGIDToMolID = {}
        for molID, molDict in enumerate(CGToAAIDMaster):
            for CGID in molDict.keys():
                CGIDToMolID[CGID] = molID
        interChainHop = []
        intraChainHop = []
        interChainTI = []
        intraChainTI = []
        interChainTITrim = []
        intraChainTITrim = []
        lambdaij = 0.3063 * elementaryCharge
        T = 290
        for chromophore in chromophoreList:
            mol1ID = CGIDToMolID[chromophore.CGIDs[0]]
            for index, neighbour in enumerate(chromophore.neighbours):
                chromophore2 = chromophoreList[neighbour[0]]
                mol2ID = CGIDToMolID[chromophore2.CGIDs[0]]
                hopRate = calculateHopRate(lambdaij, chromophore.neighboursTI[index] * elementaryCharge, chromophore.neighboursDeltaE[index] * elementaryCharge, T)
                if mol1ID == mol2ID:
                    if chromophore2.ID > chromophore.ID:
                        if hopRate != 0:
                            intraChainHop.append(hopRate)
                        intraChainTI.append(chromophore.neighboursTI[index])
                        if chromophore.neighboursTI[index] != 0:
                            intraChainTITrim.append(chromophore.neighboursTI[index])
                else:
                    if chromophore2.ID > chromophore.ID:
                        if hopRate != 0:
                            interChainHop.append(hopRate)
                        interChainTI.append(chromophore.neighboursTI[index])
                        if chromophore.neighboursTI[index] != 0:
                            interChainTITrim.append(chromophore.neighboursTI[index])

        plotHist(tempDir, intraChainHop, 'intraChainHop')
        plotHist(tempDir, interChainHop, 'interChainHop')
        plotHist(tempDir, intraChainTI, 'intraChainTI')
        plotHist(tempDir, interChainTI, 'interChainTI')
        plotHist(tempDir, intraChainTITrim, 'intraChainTITrim')
        plotHist(tempDir, interChainTITrim, 'interChainTITrim')

        print len(intraChainHop), len(interChainHop)
        plotHist(tempDir, intraChainHop, 'hopMix', xvals=interChainHop)
