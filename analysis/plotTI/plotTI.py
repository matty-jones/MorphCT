import sys
import os
import matplotlib
matplotlib.use('Agg')
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
        plt.xticks([1E0, 1E3, 1E6, 1E9, 1E12, 1E15, 1E18])
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
        plt.xticks([1E0, 1E3, 1E6, 1E9, 1E12, 1E15, 1E18])
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


def splitMolecules(inputDictionary):
    # Split the full morphology into individual molecules
    moleculeAAIDs = []
    moleculeLengths = []
    # Create a lookup table `neighbour list' for all connected atoms called {bondedAtoms}
    bondedAtoms = helperFunctions.obtainBondedList(inputDictionary['bond'])
    moleculeList = [i for i in range(len(inputDictionary['type']))]
    # Recursively add all atoms in the neighbour list to this molecule
    for molID in range(len(moleculeList)):
        moleculeList = updateMolecule(molID, moleculeList, bondedAtoms)
    # Create a dictionary of the molecule data
    moleculeData = {}
    for atomID in range(len(inputDictionary['type'])):
        if moleculeList[atomID] not in moleculeData:
            moleculeData[moleculeList[atomID]] = [atomID]
        else:
            moleculeData[moleculeList[atomID]].append(atomID)
    # Return the list of AAIDs and the lengths of the molecules
    for moleculeID in list(moleculeData.keys()):
        moleculeAAIDs.append(sorted(moleculeData[moleculeID]))
        moleculeLengths.append(len(moleculeData[moleculeID]))
    return moleculeAAIDs, moleculeLengths


def updateMolecule(atomID, moleculeList, bondedAtoms):
    # Recursively add all neighbours of atom number atomID to this molecule
    try:
        for bondedAtom in bondedAtoms[atomID]:
            # If the moleculeID of the bonded atom is larger than that of the current one,
            # update the bonded atom's ID to the current one's to put it in this molecule,
            # then iterate through all of the bonded atom's neighbours
            if moleculeList[bondedAtom] > moleculeList[atomID]:
                moleculeList[bondedAtom] = moleculeList[atomID]
                moleculeList = updateMolecule(bondedAtom, moleculeList, bondedAtoms)
            # If the moleculeID of the current atom is larger than that of the bonded one,
            # update the current atom's ID to the bonded one's to put it in this molecule,
            # then iterate through all of the current atom's neighbours
            elif moleculeList[bondedAtom] < moleculeList[atomID]:
                moleculeList[atomID] = moleculeList[bondedAtom]
                moleculeList = updateMolecule(atomID, moleculeList, bondedAtoms)
            # Else: both the current and the bonded atom are already known to be in this
            # molecule, so we don't have to do anything else.
    except KeyError:
        # This means that there are no bonded CG sites (i.e. it's a single molecule)
        pass
    return moleculeList


if __name__ == "__main__":
    tempDirs = []
    for fileName in os.listdir(os.getcwd()):
        if ('py' not in fileName) and ('pdf' not in fileName) and ('store' not in fileName):
            tempDirs.append(fileName)
    plt.figure()
    for tempDir in tempDirs:
        for fileName in os.listdir(os.getcwd() + '/' + tempDir):
            if ("pickle" in fileName) and (tempDir in fileName):
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

        print("Determining molecule IDs...")
        CGIDToMolID = {}
        if CGToAAIDMaster is not None:
            # Normal operation with a CGMorphology defined (fine-graining was performed)
            for molID, molDict in enumerate(CGToAAIDMaster):
                for CGID in list(molDict.keys()):
                    CGIDToMolID[CGID] = molID
        elif (len(parameterDict['CGSiteSpecies']) == 1) and (('AARigidBodySpecies' not in parameterDict) or (len(parameterDict['AARigidBodySpecies']) == 0)):   # The not in is a catch for the old PAH systems
            print("Small-molecule system detected, assuming each chromophore is its own molecule...")
            # When CGMorphology doesn't exist, and no rigid body species have been specified, then 
            # every chromophore is its own molecule)
            for index, chromo in enumerate(chromophoreList):
                for CGID in chromo.CGIDs:
                    CGIDToMolID[CGID] = chromo.ID
        else:
            # No CGMorphology, but not small molecules either, so determine molecules based on bonds
            print("Polymeric system detected, determining molecules based on AA bonds (slow calculation)...")
            moleculeAAIDs, moleculeLengths = splitMolecules(AAMorphologyDict)
            for index, moleculeAAIDList in enumerate(moleculeAAIDs):
                for AAID in moleculeAAIDList:
                    CGIDToMolID[AAID] = index
        donorInterChainHop = []
        donorIntraChainHop = []
        donorInterChainTI = []
        donorIntraChainTI = []
        donorInterChainTITrim = []
        donorIntraChainTITrim = []
        acceptorInterChainHop = []
        acceptorIntraChainHop = []
        acceptorInterChainTI = []
        acceptorIntraChainTI = []
        acceptorInterChainTITrim = []
        acceptorIntraChainTITrim = []
        try:
            if parameterDict['reorganisationEnergyDonor'] is not None:
                donorLambdaij = parameterDict['reorganisationEnergyDonor'] * elementaryCharge
            if parameterDict['reorganisationEnergyAcceptor'] is not None:
                acceptorLambdaij = parameterDict['reorganisationEnergyAcceptor'] * elementaryCharge
        except KeyError: # Old MorphCT fix
            print("Only one reorganisation energy found, assuming donor and continuing")
            donorLambdaij = parameterDict['reorganisationEnergy'] * elementaryCharge
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
            print("Donor intra-chain hop proportion =", len(donorIntraChainHop) / (len(donorIntraChainHop) + len(donorInterChainHop)))
            plotHist(tempDir, donorIntraChainHop, 'DonorHopMix', xvals=donorInterChainHop)
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
            print("Acceptor intra-chain hop proportion =", len(acceptorIntraChainHop) / (len(acceptorIntraChainHop) + len(acceptorInterChainHop)))
            plotHist(tempDir, acceptorIntraChainHop, 'AcceptorHopMix', xvals=acceptorInterChainHop)


        if len(HOMOLevels) > 0:
            print("DONOR HOMO LEVEL =", np.average(HOMOLevels), "+/-", np.std(HOMOLevels)/np.sqrt(len(HOMOLevels)))
            print("DONOR BANDGAP =", np.average(donorBandgap), "+/-", np.std(donorBandgap)/np.sqrt(len(donorBandgap)))
        if len(LUMOLevels) > 0:
            print("ACCEPTOR LUMO LEVEL =", np.average(LUMOLevels), "+/-", np.std(LUMOLevels)/np.sqrt(len(LUMOLevels)))
            print("ACCEPTOR BANDGAP =", np.average(acceptorBandgap), "+/-", np.std(acceptorBandgap)/np.sqrt(len(acceptorBandgap)))
        print("\n")
