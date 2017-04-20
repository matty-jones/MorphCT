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


def gaussFit(data, materialType):
    n = len(data)
    mean = np.mean(data)
    std = np.std(data)
    print(materialType, "Delta Eij stats: mean =", mean, "std =", std)
    hist, binEdges = np.histogram(data, bins=100)
    try:
        fitArgs, fitConv = curve_fit(gaussian, binEdges[:-1], hist, p0=[1, mean, std])
    except RuntimeError:
        return None, None
    return binEdges, fitArgs


########## NEED TO UPDATE THIS FUNCTION WITH THE NEW MODES BELOW
def plotHist(saveDir, yvals, mode, xvals=None, gaussBins=None, fitArgs=None):
    if mode == 'HOMO':
        plt.hist(yvals, 20, color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('HOMO Level (eV)')
        plt.xlim([-6.5, -4.0])
        fileName = 'HOMODoS.pdf'

    if mode == 'LUMO':
        plt.hist(yvals, 20, color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('LUMO Level (eV)')
        plt.xlim([-4.0, -2.5])
        fileName = 'LUMODoS.pdf'

    elif mode == 'DonorBandgap':
        plt.hist(yvals, 20, color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Donor Bandgap (eV)')
        plt.xlim([3, 9])
        fileName = 'DonorBandgap.pdf'

    elif mode == 'AcceptorBandgap':
        plt.hist(yvals, 20, color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Acceptor Bandgap (eV)')
        plt.xlim([3, 9])
        fileName = 'AcceptorBandgap.pdf'

    elif mode == 'DonorTI':
        plt.hist(yvals, 20, color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Donor Transfer Integral (eV)')
        plt.xlim([0.0, 1.2])
        fileName = 'DonorTI.pdf'

    elif mode == 'AcceptorTI':
        plt.hist(yvals, 20, color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Acceptor Transfer Integral (eV)')
        plt.xlim([0.0, 1.2])
        fileName = 'AcceptorTI.pdf'

    elif mode == 'DonorTITrimmed':
        plt.hist(yvals, 20, color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Donor Non-Zero Transfer Integral (eV)')
        plt.xlim([0.0, 1.2])
        fileName = 'DonorTITrimmed.pdf'

    elif mode == 'AcceptorTITrimmed':
        plt.hist(yvals, 20, color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Acceptor Non-Zero Transfer Integral (eV)')
        plt.xlim([0.0, 1.2])
        fileName = 'AcceptorTITrimmed.pdf'

    elif mode == 'DonorDeltaEij':
        n, bins, patches = plt.hist(yvals, 20, color = ['b'])
        if gaussBins is not None:
            gaussY = gaussian(gaussBins[:-1], *fitArgs)
            scaleFactor = max(n)/max(gaussY)
            plt.plot(gaussBins[:-1], gaussY*scaleFactor, 'ro:')
        plt.ylabel('Frequency')
        plt.xlabel('Donor Delta Eij (eV)')
        plt.xlim([-0.5, 0.5])
        fileName = 'DonorDeltaEij.pdf'

    elif mode == 'AcceptorDeltaEij':
        n, bins, patches = plt.hist(yvals, 20, color = ['b'])
        if gaussBins is not None:
            gaussY = gaussian(gaussBins[:-1], *fitArgs)
            scaleFactor = max(n)/max(gaussY)
            plt.plot(gaussBins[:-1], gaussY*scaleFactor, 'ro:')
        plt.ylabel('Frequency')
        plt.xlabel('Acceptor Delta Eij (eV)')
        plt.xlim([-0.5, 0.5])
        fileName = 'AcceptorDeltaEij.pdf'

    elif mode == 'DonorIntraChainHop':
        if len(yvals) > 0:
            plt.hist(yvals, bins = np.logspace(1, 18, 40), color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Donor Intra-Chain Hop rate (s' + r'^{-1}' + ')')
        plt.gca().set_xscale('log')
        plt.xlim([1, 1E18])
        fileName = 'DonorIntrakij.pdf'

    elif mode == 'DonorInterChainHop':
        if len(yvals) > 0:
            plt.hist(yvals, bins = np.logspace(1, 18, 40), color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Donor Inter-Chain Hop rate (s' + r'$^{-1}$' + ')')
        plt.gca().set_xscale('log')
        plt.xlim([1, 1E18])
        fileName = 'DonorInterkij.pdf'

    elif mode == 'DonorIntraChainTI':
        if len(yvals) > 0:
            plt.hist(yvals, 20, color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Donor Intra-Chain TI (eV)')
        plt.xlim([0, 1.2])
        fileName = 'DonorIntraTij.pdf'

    elif mode == 'DonorInterChainTI':
        if len(yvals) > 0:
            plt.hist(yvals, 20, color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Donor Inter-Chain TI (eV)')
        plt.xlim([0, 1.2])
        fileName = 'DonorInterTij.pdf'

    elif mode == 'DonorIntraChainTITrim':
        if len(yvals) > 0:
            plt.hist(yvals, 20, color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Donor Intra-Chain TI (eV)')
        plt.xlim([0, 1.2])
        fileName = 'DonorIntraTijTrim.pdf'

    elif mode == 'DonorInterChainTITrim':
        if len(yvals) > 0:
            plt.hist(yvals, 20, color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Donor Inter-Chain TI (eV)')
        plt.xlim([0, 1.2])
        fileName = 'DonorInterTijTrim.pdf'

    elif mode == 'DonorHopMix':
        plt.hist([yvals, xvals], bins = np.logspace(1, 18, 40), stacked = True, color = ['r', 'b'], label = ['Intra-Molecular', 'Inter-Molecular'])
        plt.ylabel('Frequency')
        plt.xlabel('Donor Hopping Rate (s' + r'$^{-1}$' + ')')
        plt.xlim([1,1E18])
        #plt.ylim([0,8000])
        plt.legend(loc = 2, prop = {'size':18})
        plt.gca().set_xscale('log')
        fileName = 'DonorHoppingRateMixed.pdf'

    elif mode == 'AcceptorIntraChainHop':
        if len(yvals) > 0:
            plt.hist(yvals, bins = np.logspace(1, 18, 40), color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Acceptor Intra-Chain Hop rate (s' + r'^{-1}' + ')')
        plt.gca().set_xscale('log')
        plt.xlim([1, 1E18])
        fileName = 'AcceptorIntrakij.pdf'

    elif mode == 'AcceptorInterChainHop':
        if len(yvals) > 0:
            plt.hist(yvals, bins = np.logspace(1, 18, 40), color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Acceptor Inter-Chain Hop rate (s' + r'$^{-1}$' + ')')
        plt.gca().set_xscale('log')
        plt.xlim([1, 1E18])
        fileName = 'AcceptorInterkij.pdf'

    elif mode == 'AcceptorIntraChainTI':
        if len(yvals) > 0:
            plt.hist(yvals, 20, color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Acceptor Intra-Chain TI (eV)')
        plt.xlim([0, 1.2])
        fileName = 'AcceptorIntraTij.pdf'

    elif mode == 'AcceptorInterChainTI':
        if len(yvals) > 0:
            plt.hist(yvals, 20, color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Acceptor Inter-Chain TI (eV)')
        plt.xlim([0, 1.2])
        fileName = 'AcceptorInterTij.pdf'

    elif mode == 'AcceptorIntraChainTITrim':
        if len(yvals) > 0:
            plt.hist(yvals, 20, color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Acceptor Intra-Chain TI (eV)')
        plt.xlim([0, 1.2])
        fileName = 'AcceptorIntraTijTrim.pdf'

    elif mode == 'AcceptorInterChainTITrim':
        if len(yvals) > 0:
            plt.hist(yvals, 20, color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Acceptor Inter-Chain TI (eV)')
        plt.xlim([0, 1.2])
        fileName = 'AcceptorInterTijTrim.pdf'

    elif mode == 'AcceptorHopMix':
        plt.hist([yvals, xvals], bins = np.logspace(1, 18, 40), stacked = True, color = ['r', 'b'], label = ['Intra-Molecular', 'Inter-Molecular'])
        plt.ylabel('Frequency')
        plt.xlabel('Acceptor Hopping Rate (s' + r'$^{-1}$' + ')')
        plt.xlim([1,1E18])
        #plt.ylim([0,8000])
        plt.legend(loc = 2, prop = {'size':18})
        plt.gca().set_xscale('log')
        fileName = 'AcceptorHoppingRateMixed.pdf'

    plt.savefig(saveDir+'/'+fileName)
    plt.clf()
    print("Figure saved as", saveDir+"/"+fileName)


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
        if ('py' not in fileName) and ('pdf' not in fileName) and ('store' not in fileName):
            tempDirs.append(fileName)
    plt.figure()
    for tempDir in tempDirs:
        for fileName in os.listdir(os.getcwd() + '/' + tempDir):
            if "pickle" in fileName:
                mainMorphologyPickleName = os.getcwd() + '/' + tempDir + '/' + fileName
        AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = helperFunctions.loadPickle(mainMorphologyPickleName)
        HOMOLevels = []
        LUMOLevels = []
        donorBandgap = []
        donorTransferIntegrals = []
        donorTransferIntegralsTrimmed = []
        donorDeltaEij = []
        acceptorBandgap = []
        acceptorTransferIntegrals = []
        acceptorTransferIntegralsTrimmed = []
        acceptorDeltaEij = []
        for chromophore in chromophoreList:
            if chromophore.species == 'Donor':
                HOMOLevels.append(chromophore.HOMO)
                donorBandgap.append(chromophore.LUMO - chromophore.HOMO)
                for neighbourLoc, transferIntegral in enumerate(chromophore.neighboursTI):
                    if (transferIntegral is not None) and (chromophore.neighboursDeltaE[neighbourLoc] is not None):
                        donorTransferIntegrals.append(transferIntegral)
                        donorDeltaEij.append(chromophore.neighboursDeltaE[neighbourLoc])
                for neighbourTI in chromophore.neighboursTI:
                    if (neighbourTI != 0) and (neighbourTI is not None):
                        donorTransferIntegralsTrimmed.append(neighbourTI)
            else:
                LUMOLevels.append(chromophore.LUMO)
                acceptorBandgap.append(chromophore.LUMO - chromophore.HOMO)
                for neighbourLoc, transferIntegral in enumerate(chromophore.neighboursTI):
                    if (transferIntegral is not None) and (chromophore.neighboursDeltaE[neighbourLoc] is not None):
                        acceptorTransferIntegrals.append(transferIntegral)
                        acceptorDeltaEij.append(chromophore.neighboursDeltaE[neighbourLoc])
                for neighbourTI in chromophore.neighboursTI:
                    if (neighbourTI != 0) and (neighbourTI is not None):
                        acceptorTransferIntegralsTrimmed.append(neighbourTI)

        donorBinEdges, donorFitArgs = gaussFit(donorDeltaEij, 'Donor')
        acceptorBinEdges, acceptorFitArgs = gaussFit(acceptorDeltaEij, 'Acceptor')

        if len(HOMOLevels) > 0:
            plotHist(tempDir, HOMOLevels, 'HOMO')
            plotHist(tempDir, donorTransferIntegrals, 'DonorTI')
            plotHist(tempDir, donorTransferIntegralsTrimmed, 'DonorTITrimmed')
            plotHist(tempDir, donorDeltaEij, 'DonorDeltaEij', gaussBins=donorBinEdges, fitArgs=donorFitArgs)
            plotHist(tempDir, donorBandgap, 'DonorBandgap')
        if len(LUMOLevels) > 0:
            plotHist(tempDir, LUMOLevels, 'LUMO')
            plotHist(tempDir, acceptorTransferIntegrals, 'AcceptorTI')
            plotHist(tempDir, acceptorTransferIntegralsTrimmed, 'AcceptorTITrimmed')
            plotHist(tempDir, acceptorDeltaEij, 'AcceptorDeltaEij', gaussBins=acceptorBinEdges, fitArgs=acceptorFitArgs)
            plotHist(tempDir, acceptorBandgap, 'AcceptorBandgap')

        CGIDToMolID = {}
        # Edit for when CGToAAIDMaster doesn't exist (every chromophore is its own molecule)
        try:
            for molID, molDict in enumerate(CGToAAIDMaster):
                for CGID in list(molDict.keys()):
                    CGIDToMolID[CGID] = molID
        except TypeError:
            for index, chromo in enumerate(chromophoreList):
                for CGID in chromo.CGIDs:
                    CGIDToMolID[CGID] = chromo.ID
        donorInterChainHop = []
        donorIntraChainHop = []
        donorInterChainTI = []
        donorIntraChainTI = []
        donorInterChainTITrim = []
        donorIntraChainTITrim = []
        donorLambdaij = parameterDict['reorganisationEnergyDonor'] * elementaryCharge
        acceptorInterChainHop = []
        acceptorIntraChainHop = []
        acceptorInterChainTI = []
        acceptorIntraChainTI = []
        acceptorInterChainTITrim = []
        acceptorIntraChainTITrim = []
        acceptorLambdaij = parameterDict['reorganisationEnergyAcceptor'] * elementaryCharge
        T = 290
        for chromophore in chromophoreList:
            mol1ID = CGIDToMolID[chromophore.CGIDs[0]]
            for index, neighbour in enumerate(chromophore.neighbours):
                if chromophore.neighboursTI[index] is None:
                    continue
                chromophore2 = chromophoreList[neighbour[0]]
                mol2ID = CGIDToMolID[chromophore2.CGIDs[0]]
                if chromophore.species == 'Donor':
                    hopRate = calculateHopRate(donorLambdaij, chromophore.neighboursTI[index] * elementaryCharge, chromophore.neighboursDeltaE[index] * elementaryCharge, T)
                elif chromophore.species == 'Acceptor':
                    hopRate = calculateHopRate(acceptorLambdaij, chromophore.neighboursTI[index] * elementaryCharge, chromophore.neighboursDeltaE[index] * elementaryCharge, T)
                if mol1ID == mol2ID:
                    if chromophore2.ID > chromophore.ID:
                        if chromophore.species == 'Donor':
                            if hopRate != 0:
                                donorIntraChainHop.append(hopRate)
                            donorIntraChainTI.append(chromophore.neighboursTI[index])
                            if chromophore.neighboursTI[index] != 0:
                                donorIntraChainTITrim.append(chromophore.neighboursTI[index])
                        elif chromophore.species == 'Acceptor':
                            if hopRate != 0:
                                acceptorIntraChainHop.append(hopRate)
                            acceptorIntraChainTI.append(chromophore.neighboursTI[index])
                            if chromophore.neighboursTI[index] != 0:
                                acceptorIntraChainTITrim.append(chromophore.neighboursTI[index])
                else:
                    if chromophore2.ID > chromophore.ID:
                        if chromophore.species == 'Donor':
                            if hopRate != 0:
                                donorInterChainHop.append(hopRate)
                            donorInterChainTI.append(chromophore.neighboursTI[index])
                            if chromophore.neighboursTI[index] != 0:
                                donorInterChainTITrim.append(chromophore.neighboursTI[index])
                        if chromophore.species == 'Acceptor':
                            if hopRate != 0:
                                acceptorInterChainHop.append(hopRate)
                            acceptorInterChainTI.append(chromophore.neighboursTI[index])
                            if chromophore.neighboursTI[index] != 0:
                                acceptorInterChainTITrim.append(chromophore.neighboursTI[index])

        if len(HOMOLevels) > 0:
            if len(donorIntraChainHop) > 0:
                plotHist(tempDir, donorIntraChainHop, 'DonorIntraChainHop')
                plotHist(tempDir, donorIntraChainTI, 'DonorIntraChainTI')
                plotHist(tempDir, donorIntraChainTITrim, 'DonorIntraChainTITrim')
            if len(donorInterChainHop) > 0:
                plotHist(tempDir, donorInterChainHop, 'DonorInterChainHop')
                plotHist(tempDir, donorInterChainTI, 'DonorInterChainTI')
                plotHist(tempDir, donorInterChainTITrim, 'DonorInterChainTITrim')
            print("Donor intra-chain hops =", len(donorIntraChainHop), "Donor inter-chain hops =", len(donorInterChainHop))
        if len(LUMOLevels) > 0:
            if len(acceptorIntraChainHop) > 0:
                plotHist(tempDir, acceptorIntraChainHop, 'AcceptorIntraChainHop')
                plotHist(tempDir, acceptorIntraChainTI, 'AcceptorIntraChainTI')
                plotHist(tempDir, acceptorIntraChainTITrim, 'AcceptorIntraChainTITrim')
            if len(acceptorInterChainHop) > 0:
                plotHist(tempDir, acceptorInterChainHop, 'AcceptorInterChainHop')
                plotHist(tempDir, acceptorInterChainTI, 'AcceptorInterChainTI')
                plotHist(tempDir, acceptorInterChainTITrim, 'AcceptorInterChainTITrim')
            print("Acceptor intra-chain hops =", len(acceptorIntraChainHop), "Acceptor inter-chain hops =", len(acceptorInterChainHop))

        plotHist(tempDir, donorIntraChainHop, 'DonorHopMix', xvals=donorInterChainHop)
        plotHist(tempDir, acceptorIntraChainHop, 'AcceptorHopMix', xvals=acceptorInterChainHop)
