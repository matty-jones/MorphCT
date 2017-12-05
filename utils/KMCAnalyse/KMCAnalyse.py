import os
import sys
import pickle
import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import scipy.stats
from scipy.sparse import lil_matrix
sys.path.append('../../code/')
sys.path.append('../code')
import helperFunctions
try:
    import mpl_toolkits.mplot3d as p3
except ImportError:
    print("Could not import 3D plotting engine, calling the plotMolecule3D function will result in an error!")
from collections import OrderedDict
import shutil
import glob
import re
import argparse
import copy

elementaryCharge = 1.60217657E-19  # C
kB = 1.3806488E-23  # m^{2} kg s^{-2} K^{-1}
hbar = 1.05457173E-34 # m^{2} kg s^{-1}
temperature = 290  # K


def loadKMCResultsPickle(directory):
    try:
        with open(directory + '/KMC/KMCResults.pickle', 'rb') as pickleFile:
            carrierData = pickle.load(pickleFile)
    except FileNotFoundError:
        print("No final KMCResults.pickle found. Creating it from incomplete parts...")
        createResultsPickle(directory)
        with open(directory + '/KMC/KMCResults.pickle', 'rb') as pickleFile:
            carrierData = pickle.load(pickleFile)
    except UnicodeDecodeError:
        with open(directory + '/KMC/KMCResults.pickle', 'rb') as pickleFile:
            carrierData = pickle.load(pickleFile, encoding='latin1')
    return carrierData


def splitCarriersByType(carrierData):
    # If only one carrier type has been given, call the carriers holes and skip the electron calculations
    listVariables = ['currentTime', 'ID', 'noHops', 'displacement', 'lifetime', 'finalPosition', 'image', 'initialPosition']
    try:
        carrierDataHoles = {'carrierHistoryMatrix': carrierData['holeHistoryMatrix'], 'seed': carrierData['seed']}
        carrierDataElectrons = {'carrierHistoryMatrix': carrierData['electronHistoryMatrix'], 'seed': carrierData['seed']}
        for listVar in listVariables:
            carrierDataHoles[listVar] = []
            carrierDataElectrons[listVar] = []
            for carrierIndex, chargeType in enumerate(carrierData['carrierType']):
                if chargeType == 'Hole':
                    carrierDataHoles[listVar].append(carrierData[listVar][carrierIndex])
                elif chargeType == 'Electron':
                    carrierDataElectrons[listVar].append(carrierData[listVar][carrierIndex])
    except:
        print("Multiple charge carriers not found, assuming donor material and holes only")
        try:
            carrierDataHoles = {'carrierHistoryMatrix': carrierData['carrierHistoryMatrix'], 'seed': carrierData['seed']}
        except KeyError:
            carrierDataHoles = {'carrierHistoryMatrix': carrierData['carrierHistoryMatrix'], 'seed': 0}
        carrierDataElectrons = None
        for listVar in listVariables:
            carrierDataHoles[listVar] = []
            for carrierIndex, carrierID in enumerate(carrierData['ID']):
                carrierDataHoles[listVar].append(carrierData[listVar][carrierIndex])
    return carrierDataHoles, carrierDataElectrons


def getCarrierData(carrierData):
    try:
        carrierHistory = carrierData['carrierHistoryMatrix']
    except:
        carrierHistory = None
    totalDataPoints = 0
    totalDataPointsAveragedOver = 0
    squaredDisps = {}
    actualTimes = {}
    carrierTypes = {}
    for carrierIndex, displacement in enumerate(carrierData['displacement']):
        if (carrierData['currentTime'][carrierIndex] > carrierData['lifetime'][carrierIndex] * 2) or (carrierData['currentTime'][carrierIndex] < carrierData['lifetime'][carrierIndex] / 2.0) or (carrierData['noHops'][carrierIndex] == 1):
            totalDataPoints += 1
            continue
        carrierKey = str(carrierData['lifetime'][carrierIndex])
        if carrierKey not in squaredDisps:
            squaredDisps[carrierKey] = [(carrierData['displacement'][carrierIndex] * 1E-10) ** 2]  # Carrier displacement is in angstroems, convert to metres
            actualTimes[carrierKey] = [carrierData['currentTime'][carrierIndex]]
        else:
            squaredDisps[carrierKey].append((carrierData['displacement'][carrierIndex] * 1E-10) ** 2)  # Carrier displacement is in angstroems, convert to metres
            actualTimes[carrierKey].append(carrierData['currentTime'][carrierIndex])
        # Also keep track of whether each carrier is a hole or an electron
        totalDataPointsAveragedOver += 1
        totalDataPoints += 1
    times = []
    MSDs = []
    timeStandardErrors = []
    MSDStandardErrors = []
    for time, disps in squaredDisps.items():
        times.append(float(time))
        timeStandardErrors.append(np.std(actualTimes[time]) / len(actualTimes[time]))
        MSDs.append(np.average(disps))
        MSDStandardErrors.append(np.std(disps) / len(disps))
    return carrierHistory, times, MSDs, timeStandardErrors, MSDStandardErrors


def plotConnections(chromophoreList, simDims, carrierHistory, directory, carrierType):
    # A complicated function that shows connections between carriers in 3D that carriers prefer to hop between.
    # Connections that are frequently used are highlighted in black, whereas rarely used connections are more white.
    # Find a good normalisation factor
    carrierHistory = carrierHistory.toarray()
    normalizeTo = np.max(carrierHistory)
    # Try to get the colour map first
    colormap = plt.cm.plasma
    minimum = np.min(carrierHistory[np.nonzero(carrierHistory)])
    maximum = np.max(carrierHistory[np.nonzero(carrierHistory)])
    plt.gcf()
    levels = np.linspace(np.log10(minimum), np.log10(maximum), 100)
    coloursForMap = plt.contourf([[0, 0], [0, 0]], levels, cmap = colormap)
    plt.clf()
    # Now for the actual plot
    fig = plt.gcf()
    ax = p3.Axes3D(fig)
    for chromo1, row in enumerate(carrierHistory):
        for chromo2, value in enumerate(row):
            if value > 0:
                coords1 = chromophoreList[chromo1].posn
                coords2 = chromophoreList[chromo2].posn
                # Only plot connections between chromophores in the same image
                plotConnection = True
                for neighbour in chromophoreList[chromo1].neighbours:
                    if neighbour[0] != chromophoreList[chromo2].ID:
                        continue
                    if neighbour[1] != [0, 0, 0]:
                        plotConnection = False
                        break
                if plotConnection is True:
                    #ax.scatter(coords1[0], coords1[1], coords1[2], c = 'k', s = '5')
                    #ax.scatter(coords2[0], coords2[1], coords2[2], c = 'k', s = '5')
                    line = [coords2[0] - coords1[0], coords2[1] - coords1[1], coords2[2] - coords2[1]]
                    if (np.abs(coords2[0] - coords1[0]) < simDims[0][1]) and (np.abs(coords2[1] - coords1[1]) < simDims[1][1]) and (np.abs(coords2[2] - coords1[2]) < simDims[2][1]):
                        #colourIntensity = value / normalizeTo
                        colourIntensity = np.log10(value) / np.log10(normalizeTo)
                        ax.plot([coords1[0], coords2[0]], [coords1[1], coords2[1]], [coords1[2], coords2[2]], c = colormap(colourIntensity), linewidth = 0.5, alpha = colourIntensity)
    # Draw boxlines
    # Varying X
    ax.plot([simDims[0][0], simDims[0][1]], [simDims[1][0], simDims[1][0]], [simDims[2][0], simDims[2][0]], c = 'k', linewidth = 1.0)
    ax.plot([simDims[0][0], simDims[0][1]], [simDims[1][1], simDims[1][1]], [simDims[2][0], simDims[2][0]], c = 'k', linewidth = 1.0)
    ax.plot([simDims[0][0], simDims[0][1]], [simDims[1][0], simDims[1][0]], [simDims[2][1], simDims[2][1]], c = 'k', linewidth = 1.0)
    ax.plot([simDims[0][0], simDims[0][1]], [simDims[1][1], simDims[1][1]], [simDims[2][1], simDims[2][1]], c = 'k', linewidth = 1.0)
    # Varying Y
    ax.plot([simDims[0][0], simDims[0][0]], [simDims[1][0], simDims[1][1]], [simDims[2][0], simDims[2][0]], c = 'k', linewidth = 1.0)
    ax.plot([simDims[0][1], simDims[0][1]], [simDims[1][0], simDims[1][1]], [simDims[2][0], simDims[2][0]], c = 'k', linewidth = 1.0)
    ax.plot([simDims[0][0], simDims[0][0]], [simDims[1][0], simDims[1][1]], [simDims[2][1], simDims[2][1]], c = 'k', linewidth = 1.0)
    ax.plot([simDims[0][1], simDims[0][1]], [simDims[1][0], simDims[1][1]], [simDims[2][1], simDims[2][1]], c = 'k', linewidth = 1.0)
    # Varying Z
    ax.plot([simDims[0][0], simDims[0][0]], [simDims[1][0], simDims[1][0]], [simDims[2][0], simDims[2][1]], c = 'k', linewidth = 1.0)
    ax.plot([simDims[0][0], simDims[0][0]], [simDims[1][1], simDims[1][1]], [simDims[2][0], simDims[2][1]], c = 'k', linewidth = 1.0)
    ax.plot([simDims[0][1], simDims[0][1]], [simDims[1][0], simDims[1][0]], [simDims[2][0], simDims[2][1]], c = 'k', linewidth = 1.0)
    ax.plot([simDims[0][1], simDims[0][1]], [simDims[1][1], simDims[1][1]], [simDims[2][0], simDims[2][1]], c = 'k', linewidth = 1.0)

    tickLocation = range(0, int(np.log10(maximum)) + 1, 1)
    cbar = plt.colorbar(coloursForMap, ticks=tickLocation)#np.linspace(np.log10(minimum), np.log10(maximum), 6))
    cbar.ax.set_yticklabels([r'10$^{{{}}}$'.format(x) for x in tickLocation])
    plt.title('Network (' + carrierType + ')', y = 1.1)
    fileName = '01_3d' + carrierType + '.pdf'
    plt.savefig(directory + '/figures/' + fileName, bbox_inches='tight')
    print("Figure saved as", directory + "/figures/" + fileName)
    plt.clf()


def calcMobility(linFitX, linFitY, avTimeError, avMSDError):
    # YVals have a std error avMSDError associated with them
    # XVals have a std error avTimeError assosciated with them
    numerator = linFitY[-1] - linFitY[0]
    denominator = linFitX[-1] - linFitX[0]
    diffusionCoeff = numerator / denominator
    # The error in the mobility is the proportionally the same as the error in the diffusion coefficient as the other variables are constants with zero error
    diffError = diffusionCoeff * np.sqrt((avMSDError / numerator)**2 + (avTimeError / denominator)**2)
    # Use Einstein relation (include the factor of 1/6!! It is in the Carbone/Troisi 2014 paper)
    mobility = elementaryCharge*diffusionCoeff/(6*kB*temperature) # This is in m^{2} / Vs
    # Convert to cm^{2}/ Vs
    mobility *= (100**2)
    mobError = (diffError / diffusionCoeff) * mobility
    return mobility, mobError


def plotMSD(times, MSDs, timeStandardErrors, MSDStandardErrors, directory, carrierType):
    ### DEBUG TEST ###
    #print "DEBUG TEST CODE ACTIVE, DELETE TO GET PROPER RESULTS!"
    #times = times[-3:]
    #MSDs = MSDs[-3:]
    #timeStandardErrors = timeStandardErrors[-3:]
    #MSDStandardErrors = MSDStandardErrors[-3:]
    ##################
    fit = np.polyfit(times, MSDs, 1)
    fitX = np.linspace(np.min(times), np.max(times), 100)
    gradient, intercept, rVal, pVal, stdErr = scipy.stats.linregress(times, MSDs)
    print("StandardError", stdErr)
    print("Fitting rVal =", rVal)
    fitY = (fitX * gradient) + intercept
    mobility, mobError = calcMobility(fitX, fitY, np.average(timeStandardErrors), np.average(MSDStandardErrors))
    plt.plot(times, MSDs)
    plt.errorbar(times, MSDs, xerr = timeStandardErrors, yerr = MSDStandardErrors)
    plt.plot(fitX, fitY, 'r')
    plt.xlabel('Time (s)')
    plt.ylabel('MSD (m'+r'$^{2}$)')
    mobilityString = '%.3e' % mobility
    #plt.title(r'$\mu_{0}$' + ' ' + carrierType + ' = ' + mobilityString + ' cm' + r'$^{2}$/Vs' % (mobility), y = 1.1)
    plt.title(r'$\mu_{0,' + carrierType[0] + r'}$' + ' = ' + mobilityString + ' cm' + r'$^{2}$/Vs' % (mobility), y = 1.1)
    fileName = '18_LinMSD' + carrierType + '.pdf'
    plt.savefig(directory + '/figures/' + fileName)
    plt.clf()
    print("Figure saved as", directory + "/figures/" + fileName)
    plt.semilogx(times, MSDs)
    plt.errorbar(times, MSDs, xerr = timeStandardErrors, yerr = MSDStandardErrors)
    plt.semilogx(fitX, fitY, 'r')
    plt.xlabel('Time (s)')
    plt.ylabel('MSD (m'+r'$^{2}$)')
    mobilityString = '%.3e' % mobility
    plt.title(r'$\mu_{0,' + carrierType[0] + r'}$' + ' = ' + mobilityString + ' cm' + r'$^{2}$/Vs' % (mobility), y = 1.1)
    fileName = '19_SemiLogMSD' + carrierType + '.pdf'
    plt.savefig(directory + '/figures/' + fileName)
    plt.clf()
    print("Figure saved as", directory + "/figures/" + fileName)
    plt.plot(times, MSDs)
    plt.errorbar(times, MSDs, xerr = timeStandardErrors, yerr = MSDStandardErrors)
    plt.plot(fitX, fitY, 'r')
    plt.xlabel('Time (s)')
    plt.ylabel('MSD (m'+r'$^{2}$)')
    plt.xscale('log')
    plt.yscale('log')
    mobilityString = '%.3e' % mobility
    plt.title(r'$\mu_{0,' + carrierType[0] + r'}$' + ' = ' + mobilityString + ' cm' + r'$^{2}$/Vs' % (mobility), y = 1.1)
    fileName = '20_LogMSD' + carrierType + '.pdf'
    plt.savefig(directory + '/figures/' + fileName)
    plt.clf()
    print("Figure saved as", directory + "/figures/" + fileName)
    return mobility, mobError, rVal**2


def calculateAnisotropy(xvals, yvals, zvals):
    # First calculate the `centre of position' for the particles
    centre = [np.mean(xvals), np.mean(yvals), np.mean(zvals)]
    # First calculate the gyration tensor:
    Sxx = 0
    Sxy = 0
    Sxz = 0
    Syy = 0
    Syz = 0
    Szz = 0
    for carrierID, rawXval in enumerate(xvals):
        xval = rawXval - centre[0]
        yval = yvals[carrierID] - centre[1]
        zval = zvals[carrierID] - centre[2]
        Sxx += xval * xval
        Sxy += xval * yval
        Sxz += xval * zval
        Syy += yval * yval
        Syz += yval * zval
        Szz += zval * zval
    S = np.array([[Sxx, Sxy, Sxz], [Sxy, Syy, Syz], [Sxz, Syz, Szz]])
    eigenValues, eigenVectors = np.linalg.eig(S)
    # Diagonalisation of S is the diagonal matrix of the eigenvalues in ascending order
    # diagonalMatrix = np.diag(sorted(eigenValues))
    # We only need the eigenvalues though, no more matrix multiplication
    diagonal = sorted(eigenValues)
    # Then calculate the relative shape anisotropy (kappa**2)
    anisotropy = (3/2) * (((diagonal[0] ** 2) + (diagonal[1] ** 2) + (diagonal[2] ** 2)) / ((diagonal[0] + diagonal[1] + diagonal[2]) ** 2)) - (1/2)
    return anisotropy


def plotAnisotropy(carrierData, directory, simDims, carrierType, plot3DGraphs):
    simExtent = [value[1] - value[0] for value in simDims]
    xvals = []
    yvals = []
    zvals = []
    colours = []
    simDimsnm = list(map(list, np.array(simDims) / 10.))
    # Get the indices of the carriers that travelled the furthest
    if len(carrierData['finalPosition']) <= 1000:
        carrierIndicesToUse = range(carrierData['finalPosition'])
    else:
        displacements = copy.deepcopy(np.array(carrierData['displacement']))
        carrierIndicesToUse = displacements.argsort()[-1000:][::-1]
    for carrierNo in carrierIndicesToUse:
        posn = carrierData['finalPosition'][carrierNo]
        #if bool(sum([x < -3 or x > 3 for x in image])):
        #    continue
        position = [0.0, 0.0, 0.0]
        for axis in range(len(posn)):
            position[axis] = (carrierData['image'][carrierNo][axis] * simExtent[axis]) + posn[axis]
        xvals.append(position[0]/10.)
        yvals.append(position[1]/10.)
        zvals.append(position[2]/10.)
        colours.append('b')
    anisotropy = calculateAnisotropy(xvals, yvals, zvals)
    if not plot3DGraphs:
        return anisotropy
    print("----------====================----------")
    print(carrierType + " charge transport anisotropy calculated as", anisotropy)
    print("----------====================----------")
    # Reduce number of plot markers
    fig = plt.gcf()
    ax = p3.Axes3D(fig)
    if len(xvals) > 1000:
        xvals = xvals[0:len(xvals):len(xvals)//1000]
        yvals = yvals[0:len(yvals):len(yvals)//1000]
        zvals = zvals[0:len(zvals):len(zvals)//1000]
    plt.scatter(xvals, yvals, zs = zvals, c = colours, s = 20)
    plt.scatter(0, 0, zs = 0, c = 'r', s = 50)
    # Draw boxlines
    # Varying X
    ax.plot([simDimsnm[0][0], simDimsnm[0][1]], [simDimsnm[1][0], simDimsnm[1][0]], [simDimsnm[2][0], simDimsnm[2][0]], c = 'k', linewidth = 1.0)
    ax.plot([simDimsnm[0][0], simDimsnm[0][1]], [simDimsnm[1][1], simDimsnm[1][1]], [simDimsnm[2][0], simDimsnm[2][0]], c = 'k', linewidth = 1.0)
    ax.plot([simDimsnm[0][0], simDimsnm[0][1]], [simDimsnm[1][0], simDimsnm[1][0]], [simDimsnm[2][1], simDimsnm[2][1]], c = 'k', linewidth = 1.0)
    ax.plot([simDimsnm[0][0], simDimsnm[0][1]], [simDimsnm[1][1], simDimsnm[1][1]], [simDimsnm[2][1], simDimsnm[2][1]], c = 'k', linewidth = 1.0)
    # Varying Y
    ax.plot([simDimsnm[0][0], simDimsnm[0][0]], [simDimsnm[1][0], simDimsnm[1][1]], [simDimsnm[2][0], simDimsnm[2][0]], c = 'k', linewidth = 1.0)
    ax.plot([simDimsnm[0][1], simDimsnm[0][1]], [simDimsnm[1][0], simDimsnm[1][1]], [simDimsnm[2][0], simDimsnm[2][0]], c = 'k', linewidth = 1.0)
    ax.plot([simDimsnm[0][0], simDimsnm[0][0]], [simDimsnm[1][0], simDimsnm[1][1]], [simDimsnm[2][1], simDimsnm[2][1]], c = 'k', linewidth = 1.0)
    ax.plot([simDimsnm[0][1], simDimsnm[0][1]], [simDimsnm[1][0], simDimsnm[1][1]], [simDimsnm[2][1], simDimsnm[2][1]], c = 'k', linewidth = 1.0)
    # Varying Z
    ax.plot([simDimsnm[0][0], simDimsnm[0][0]], [simDimsnm[1][0], simDimsnm[1][0]], [simDimsnm[2][0], simDimsnm[2][1]], c = 'k', linewidth = 1.0)
    ax.plot([simDimsnm[0][0], simDimsnm[0][0]], [simDimsnm[1][1], simDimsnm[1][1]], [simDimsnm[2][0], simDimsnm[2][1]], c = 'k', linewidth = 1.0)
    ax.plot([simDimsnm[0][1], simDimsnm[0][1]], [simDimsnm[1][0], simDimsnm[1][0]], [simDimsnm[2][0], simDimsnm[2][1]], c = 'k', linewidth = 1.0)
    ax.plot([simDimsnm[0][1], simDimsnm[0][1]], [simDimsnm[1][1], simDimsnm[1][1]], [simDimsnm[2][0], simDimsnm[2][1]], c = 'k', linewidth = 1.0)
    ax.set_xlabel('X (nm)', fontsize = 20, labelpad = 40)
    ax.set_ylabel('Y (nm)', fontsize = 20, labelpad = 40)
    ax.set_zlabel('Z (nm)', fontsize = 20, labelpad = 40)
    maximum = max([max(xvals), max(yvals), max(zvals)])
    ax.set_xlim([-maximum, maximum])
    ax.set_ylim([-maximum, maximum])
    ax.set_zlim([-maximum, maximum])
    for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks() + ax.zaxis.get_major_ticks():
        tick.label.set_fontsize(16)
    #try:
    #    plt.title(carrierType + ' transport for:' + directory[directory.index('T'):directory.index('T')+directory[directory.index('T'):].index('-')], fontsize = 24)
    #except:
    #    plt.title(carrierType + ' transport for:' + directory, fontsize = 24)
    ax.dist = 11
    if carrierType == 'Hole':
        figureIndex = '08'
    elif carrierType == 'Electron':
        figureIndex = '09'
    plt.title('Anisotropy (' + carrierType + ')', y = 1.1)
    fileName = directory + '/figures/' + figureIndex + '_anisotropy' + carrierType + '.pdf'
    plt.savefig(fileName, bbox_inches='tight')
    plt.clf()
    print("Figure saved as", fileName)
    return anisotropy


def getTempVal(string):
    hyphenList = helperFunctions.findIndex(string, '-')
    tempVal = float(string[hyphenList[-2] + 2 : hyphenList[-1]])
    return tempVal


def getFrameVal(string):
    hyphenList = helperFunctions.findIndex(string, '-')
    tempVal = int(string[hyphenList[0]+1:hyphenList[1]])
    return tempVal


def plotTemperatureProgression(tempData, mobilityData, anisotropyData, carrierType, xLabel):
    plt.gcf()
    xvals = tempData
    # DEBUG
    #xvals[-1] = 1000
    yvals = list(np.array(mobilityData)[:,0])
    yerrs = list(np.array(mobilityData)[:,1])
    plt.xlabel(xLabel)
    plt.ylabel('Mobility, cm'+r'$^{2}$ '+'V'+r'$^{-1}$'+r's$^{-1}$')
    plt.title('p1-L15-f0.0-P0.1-TX.X-e0.1', fontsize = 24)
    #plt.xlim([1.4, 2.6])
    plt.semilogy(xvals, yvals, c = 'b')
    #plt.gca().set_xscale('log')
    plt.errorbar(xvals, yvals, xerr = 0, yerr = yerrs)
    fileName = './mobility' + carrierType + '.pdf'
    plt.savefig(fileName)
    plt.clf()
    print("Figure saved as " + fileName)

    plt.plot(tempData, anisotropyData, c = 'r')
    fileName = './anisotropy' + carrierType + '.pdf'
    plt.xlabel(xLabel)
    plt.ylabel(r'$\kappa$'+', Arb. U')
    plt.savefig(fileName)
    plt.clf()
    print("Figure saved as " + fileName)


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
    hist, binEdges = np.histogram(data, bins=100)
    try:
        fitArgs, fitConv = curve_fit(gaussian, binEdges[:-1], hist, p0=[1, mean, std])
    except RuntimeError:
        fitArgs = None
    return binEdges, fitArgs, mean, std


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


def getNeighbourCutOff(chromophoreList, morphologyShape, outputDir, periodic=True):
    separationDistDonor = []
    separationDistAcceptor = []
    for chromo1 in chromophoreList:
        for chromo2Details in chromo1.neighbours:
            if (chromo2Details is None) or ((periodic is False) and (not np.array_equal(chromo2Details[1], [0, 0, 0]))) or (chromo1.ID == chromophoreList[chromo2Details[0]].ID):
                continue
            chromo2 = chromophoreList[chromo2Details[0]]
            separation = np.linalg.norm((np.array(chromo2.posn) + (np.array(chromo2Details[1]) * np.array(morphologyShape))) - chromo1.posn)
            if chromo1.species == 'Donor':
                separationDistDonor.append(separation)
            elif chromo1.species == 'Acceptor':
                separationDistAcceptor.append(separation)
    cutOffs = []
    material = ['Donor', 'Acceptor']
    for materialType, separationDist in enumerate([separationDistDonor, separationDistAcceptor]):
        if len(separationDist) == 0:
            cutOffs.append(None)
            continue
        plt.figure()
        (n, binEdges, patches) = plt.hist(separationDist, bins = 20, color = 'b')
        plt.xlabel(material[materialType] + " Chromophore Separation (Ang)")
        plt.ylabel("Frequency (Arb. U.)")
        plt.savefig(outputDir + "/03_neighbourHist" + material[materialType] + ".pdf")
        plt.close()
        print("Neighbour histogram figure saved as", outputDir + "/03_neighbourHist" + material[materialType] + ".pdf")
        bins = 0.5*(binEdges[1:]+binEdges[:-1])
        bins = np.insert(bins, 0, 0)
        n = np.insert(n, 0, 0)
        dn = np.diff(n)
        minimaIndices = []
        maximaIndices = []
        previousValue = 1E99
        for index, val in enumerate(dn):
            if (previousValue <= 0) and (val > 0):
                minimaIndices.append(index)
            if (previousValue >= 0) and (val < 0):
                maximaIndices.append(index)
            previousValue = val
        # Minimum is half way between the first maximum and the first minimum of the distribution
        cutOffs.append((bins[maximaIndices[0]] + bins[minimaIndices[0]]) / 2.0)
    return cutOffs


def getStacks(chromophoreList, morphologyShape, cutOffDonor, cutOffAcceptor, periodic=True):
    cutOffs = [cutOffDonor, cutOffAcceptor]
    materialsToCheck = ['Donor', 'Acceptor']
    stackDicts = []
    for typeIndex, materialType in enumerate(materialsToCheck):
        cutOff = cutOffs[typeIndex]
        if cutOff is None:
            stackDicts.append(None)
            continue
        chromoIDs = [chromo.ID for chromo in chromophoreList if chromo.species == materialType]
        # Create a neighbourlist based on the cutoff
        neighbourDict = createNeighbourList(chromophoreList, morphologyShape, cutOff, periodic, materialType)
        # Do the usual stackList neighbourList stuff
        stackList = [_ for _ in range(len(chromophoreList))]
        for stackID in range(len(stackList)):
            stackList = updateStack(stackID, stackList, neighbourDict)
        actualStacks = [stackList[x] for x in chromoIDs]
        print("There are", len(set(actualStacks)), materialType, "stacks in the system")
        stackDict = {}
        for index, chromophore in enumerate(chromophoreList):
            if chromophore.species != materialType:
                continue
            stackDict[chromophore.ID] = stackList[index]
        stackDicts.append(stackDict)
    return stackDicts


def createNeighbourList(chromophoreList, morphologyShape, cutOff, periodic, materialType):
    neighbourDict = {}
    for chromo1 in chromophoreList:
        for [chromo2ID, relImage] in chromo1.neighbours:
            if periodic is False:
                if not np.array_equal(relImage, [0, 0, 0]):
                    continue
            if chromo1.species != materialType:
                continue
            chromo1Posn = chromo1.posn
            chromo2Posn = np.array(chromophoreList[chromo2ID].posn) + (np.array(relImage) * np.array(morphologyShape))
            separation = np.linalg.norm(chromo2Posn - chromo1Posn)
            if separation < cutOff:
                if chromo1.ID in neighbourDict.keys():
                    neighbourDict[chromo1.ID].append(chromo2ID)
                else:
                    neighbourDict[chromo1.ID] = [chromo2ID]
    return neighbourDict


def updateStack(atomID, clusterList, neighbourDict):
    try:
        for neighbour in neighbourDict[atomID]:
            if clusterList[neighbour] > clusterList[atomID]:
                clusterList[neighbour] = clusterList[atomID]
                clusterList = updateStack(neighbour, clusterList, neighbourDict)
            elif clusterList[neighbour] < clusterList[atomID]:
                clusterList[atomID] = clusterList[neighbour]
                clusterList = updateStack(neighbour, clusterList, neighbourDict)
    except KeyError:
        pass
    return clusterList


def plotStacks3D(outputDir, chromophoreList, stackDicts, simDims):
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    colours = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'w']
    stackDict = {}
    for dictionary in stackDicts:
        if dictionary is not None:
            stackDict.update(dictionary)
    stackList = {}
    for chromophore in chromophoreList:
        stackID = stackDict[chromophore.ID]
        if stackID not in stackList.keys():
            stackList[stackID] = [chromophore]
        else:
            stackList[stackID].append(chromophore)
    for stackID, chromos in enumerate(stackList.values()):
        for chromo in chromos:
            if chromo.species == 'Donor':
                ax.scatter(chromo.posn[0], chromo.posn[1], chromo.posn[2], facecolors = 'w', edgecolors = colours[stackID%8], s = 40)
            elif chromo.species == 'Acceptor':
                ax.scatter(chromo.posn[0], chromo.posn[1], chromo.posn[2], c = colours[stackID%8], edgecolors = None, s = 40)
    # Draw boxlines
    # Varying X
    ax.plot([simDims[0][0], simDims[0][1]], [simDims[1][0], simDims[1][0]], [simDims[2][0], simDims[2][0]], c = 'k', linewidth = 1.0)
    ax.plot([simDims[0][0], simDims[0][1]], [simDims[1][1], simDims[1][1]], [simDims[2][0], simDims[2][0]], c = 'k', linewidth = 1.0)
    ax.plot([simDims[0][0], simDims[0][1]], [simDims[1][0], simDims[1][0]], [simDims[2][1], simDims[2][1]], c = 'k', linewidth = 1.0)
    ax.plot([simDims[0][0], simDims[0][1]], [simDims[1][1], simDims[1][1]], [simDims[2][1], simDims[2][1]], c = 'k', linewidth = 1.0)
    # Varying Y
    ax.plot([simDims[0][0], simDims[0][0]], [simDims[1][0], simDims[1][1]], [simDims[2][0], simDims[2][0]], c = 'k', linewidth = 1.0)
    ax.plot([simDims[0][1], simDims[0][1]], [simDims[1][0], simDims[1][1]], [simDims[2][0], simDims[2][0]], c = 'k', linewidth = 1.0)
    ax.plot([simDims[0][0], simDims[0][0]], [simDims[1][0], simDims[1][1]], [simDims[2][1], simDims[2][1]], c = 'k', linewidth = 1.0)
    ax.plot([simDims[0][1], simDims[0][1]], [simDims[1][0], simDims[1][1]], [simDims[2][1], simDims[2][1]], c = 'k', linewidth = 1.0)
    # Varying Z
    ax.plot([simDims[0][0], simDims[0][0]], [simDims[1][0], simDims[1][0]], [simDims[2][0], simDims[2][1]], c = 'k', linewidth = 1.0)
    ax.plot([simDims[0][0], simDims[0][0]], [simDims[1][1], simDims[1][1]], [simDims[2][0], simDims[2][1]], c = 'k', linewidth = 1.0)
    ax.plot([simDims[0][1], simDims[0][1]], [simDims[1][0], simDims[1][0]], [simDims[2][0], simDims[2][1]], c = 'k', linewidth = 1.0)
    ax.plot([simDims[0][1], simDims[0][1]], [simDims[1][1], simDims[1][1]], [simDims[2][0], simDims[2][1]], c = 'k', linewidth = 1.0)
    ax.set_xlim([simDims[0][0], simDims[0][1]])
    ax.set_ylim([simDims[1][0], simDims[1][1]])
    ax.set_zlim([simDims[2][0], simDims[2][1]])
    plt.savefig(outputDir + "/02_stacks.pdf", bbox_inches='tight')
    plt.close()
    print("3D Stack figure saved as", outputDir + "/02_stacks.pdf")
    #plt.show()


def determineMoleculeIDs(CGToAAIDMaster, AAMorphologyDict, parameterDict, chromophoreList):
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
    return CGIDToMolID


def plotEnergyLevels(outputDir, chromophoreList, dataDict):
    HOMOLevels = []
    LUMOLevels = []
    donorDeltaEij = []
    acceptorDeltaEij = []
    for chromo in chromophoreList:
        if chromo.species == 'Donor':
            HOMOLevels.append(chromo.HOMO)
            for neighbourIndex, deltaEij in enumerate(chromo.neighboursDeltaE):
                if (deltaEij is not None) and (chromo.neighboursTI[neighbourIndex] is not None):
                    donorDeltaEij.append(deltaEij)
        else:
            LUMOLevels.append(chromo.LUMO)
            for neighbourIndex, deltaEij in enumerate(chromo.neighboursDeltaE):
                if (deltaEij is not None) and (chromo.neighboursTI[neighbourIndex] is not None):
                    acceptorDeltaEij.append(deltaEij)
    if len(donorDeltaEij) > 0:
        donorBinEdges, donorFitArgs, donorMean, donorSTD = gaussFit(donorDeltaEij)
        dataDict['donor_delta_Eij_mean'] = donorMean
        dataDict['donor_delta_Eij_std'] = donorSTD
        dataDict['donor_delta_Eij_err'] = donorSTD / np.sqrt(len(donorDeltaEij))
        HOMOAv = np.average(HOMOLevels)
        HOMOStd = np.std(HOMOLevels)
        HOMOErr = HOMOStd / np.sqrt(len(HOMOLevels))
        dataDict['donor_frontierMO_mean'] = HOMOAv
        dataDict['donor_frontierMO_std'] = HOMOStd
        dataDict['donor_frontierMO_err'] = HOMOErr
        print("Donor HOMO Level =", HOMOAv, "+/-", HOMOErr)
        print("Donor Delta Eij stats: mean =", donorMean, "+/-", donorSTD / np.sqrt(len(donorDeltaEij)))
        plotDeltaEij(donorDeltaEij, donorBinEdges, donorFitArgs, 'Donor', outputDir + '/05_DonorDeltaEij.pdf')
    if len(acceptorDeltaEij) > 0:
        acceptorBinEdges, acceptorFitArgs, acceptorMean, acceptorSTD = gaussFit(acceptorDeltaEij)
        dataDict['acceptor_delta_Eij_mean'] = acceptorMean
        dataDict['acceptor_delta_Eij_std'] = acceptorSTD
        dataDict['acceptor_delta_Eij_err'] = acceptorSTD / np.sqrt(len(acceptorDeltaEij))
        LUMOAv = np.average(LUMOLevels)
        LUMOStd = np.std(LUMOLevels)
        LUMOErr = LUMOStd / np.sqrt(len(LUMOLevels))
        dataDict['acceptor_frontierMO_mean'] = LUMOAv
        dataDict['acceptor_frontierMO_std'] = LUMOStd
        dataDict['acceptor_frontierMO_err'] = LUMOErr
        print("Acceptor LUMO Level =", LUMOAv, "+/-", LUMOErr)
        print("Acceptor Delta Eij stats: mean =", acceptorMean, "+/-", acceptorSTD / np.sqrt(len(acceptorDeltaEij)))
        plotDeltaEij(acceptorDeltaEij, acceptorBinEdges, acceptorFitArgs, 'Acceptor', outputDir + '/06_AcceptorDeltaEij.pdf')
    return dataDict


def generateDataDict():
        materials = ['donor', 'acceptor']
        materialInspecificProperties = ['name', 'density']
        noErrorProperties = ['anisotropy', 'mobility', 'mobility_rSquared']
        hopTypes = ['intra', 'inter']
        hopTargets = ['mol', 'stack']
        hopDependentProperties = ['hops', 'proportion']
        errorProperties = ['frontierMO', 'deltaEij']
        dictionaryElements = [(prop, '---') for prop in materialInspecificProperties]
        dictionaryElements += [(material + '_' + noErrorProperty, '---') for material in materials for noErrorProperty in noErrorProperties]
        dictionaryElements += [(material + '_' + hopType + '_' + hopTarget + '_' + hopProperty, '---') for material in materials for hopType in hopTypes for hopTarget in hopTargets for hopProperty in hopDependentProperties]
        dictionaryElements += [(material + '_' + errorProperty + '_' + stat, '---') for material in materials for errorProperty in errorProperties for stat in ['mean', 'std', 'err']]
        dataDict = OrderedDict(dictionaryElements)
        return dataDict


def plotDeltaEij(deltaEij, gaussBins, fitArgs, dataType, fileName):
    plt.figure()
    n, bins, patches = plt.hist(deltaEij, np.linspace(-0.5,0.5,20), color = ['b'])
    if fitArgs is not None:
        gaussY = gaussian(gaussBins[:-1], *fitArgs)
        scaleFactor = max(n)/max(gaussY)
        plt.plot(gaussBins[:-1], gaussY*scaleFactor, 'ro:')
    else:
        print("No Gaussian found (probably zero-width delta function")
    plt.ylabel('Frequency (Arb. U.)')
    plt.xlabel(dataType + ' Delta Eij (eV)')
    plt.xlim([-0.5, 0.5])
    plt.savefig(fileName)
    plt.close()
    print("Figure saved as", fileName)


def plotMixedHoppingRates(outputDir, chromophoreList, parameterDict, stackDicts, CGToMolID, dataDict):
    # Create all the empty lists we need
    hopTypes = ['intra', 'inter']
    hopTargets = ['Stack', 'Mol']
    hopProperties = ['Rates', 'TIs']
    chromoSpecies = ['Donor', 'Acceptor']
    propertyLists = {}
    for propertyName in [hopType + hopTarget + hopProperty + species for hopType in hopTypes for hopTarget in hopTargets for hopProperty in hopProperties for species in chromoSpecies]:
        propertyLists[propertyName] = []
    try:
        if parameterDict['reorganisationEnergyDonor'] is not None:
            donorLambdaij = parameterDict['reorganisationEnergyDonor']
        if parameterDict['reorganisationEnergyAcceptor'] is not None:
            acceptorLambdaij = parameterDict['reorganisationEnergyAcceptor']
    except KeyError: # Old MorphCT fix
        print("Only one reorganisation energy found, assuming donor and continuing")
        donorLambdaij = parameterDict['reorganisationEnergy']
    T = 290
    for chromo in chromophoreList:
        mol1ID = CGToMolID[chromo.CGIDs[0]]
        for index, Tij in enumerate(chromo.neighboursTI):
            if (Tij == None) or (Tij == 0):
                continue
            chromo2 = chromophoreList[chromo.neighbours[index][0]]
            mol2ID = CGToMolID[chromo2.CGIDs[0]]
            deltaE = chromo.neighboursDeltaE[index]
            if chromo.species == 'Acceptor':
                rate = calculateHopRate(acceptorLambdaij * elementaryCharge, Tij * elementaryCharge, deltaE * elementaryCharge, T)
            else:
                rate = calculateHopRate(donorLambdaij * elementaryCharge, Tij * elementaryCharge, deltaE * elementaryCharge, T)
            #try:
            if chromo2.ID < chromo.ID:
                continue
            # Do intra- / inter- stacks
            if chromo.species == 'Acceptor':
                if stackDicts[1][chromo.ID] == stackDicts[1][chromo.neighbours[index][0]]:
                    propertyLists['intraStackRatesAcceptor'].append(rate)
                    propertyLists['intraStackTIsAcceptor'].append(Tij)
                else:
                    propertyLists['interStackRatesAcceptor'].append(rate)
                    propertyLists['interStackTIsAcceptor'].append(Tij)
            else:
                if stackDicts[0][chromo.ID] == stackDicts[0][chromo.neighbours[index][0]]:
                    propertyLists['intraStackRatesDonor'].append(rate)
                    propertyLists['intraStackTIsDonor'].append(Tij)
                else:
                    propertyLists['interStackRatesDonor'].append(rate)
                    propertyLists['interStackTIsDonor'].append(Tij)
            # Now do intra- / inter- molecules
            if mol1ID == mol2ID:
                if chromo.species == 'Acceptor':
                    propertyLists['intraMolRatesAcceptor'].append(rate)
                    propertyLists['intraMolTIsAcceptor'].append(Tij)
                else:
                    propertyLists['intraMolRatesDonor'].append(rate)
                    propertyLists['intraMolTIsDonor'].append(Tij)
            else:
                if chromo.species == 'Acceptor':
                    propertyLists['interMolRatesAcceptor'].append(rate)
                    propertyLists['interMolTIsAcceptor'].append(Tij)
                else:
                    propertyLists['interMolRatesDonor'].append(rate)
                    propertyLists['interMolTIsDonor'].append(Tij)
            #except TypeError:
            #    print(repr(sys.exc_info()))
            #    print("TYPE ERROR EXCEPTION")
            #    pass
    #print(len(propertyLists['intraStackRatesDonor']), len(propertyLists['intraStackRatesAcceptor']), len(propertyLists['intraMolRatesDonor']), len(propertyLists['intraMolRatesAcceptor']))
    # Donor Stack Plots:
    if (len(propertyLists['intraStackRatesDonor']) > 0) or (len(propertyLists['interStackRatesDonor']) > 0):
        plotStackedHistRates(propertyLists['intraStackRatesDonor'], propertyLists['interStackRatesDonor'], ['Intra-Stack', 'Inter-Stack'], 'Donor', outputDir + '/16_DonorHoppingRate_Stacks.pdf')
        plotStackedHistTIs(propertyLists['intraStackTIsDonor'], propertyLists['interStackTIsDonor'], ['Intra-Stack', 'Inter-Stack'], 'Donor', outputDir + '/12_DonorTransferIntegral_Stacks.pdf')
    # Acceptor Stack Plots:
    if (len(propertyLists['intraStackRatesAcceptor']) > 0) or (len(properyLists['interStackRatesAcceptor']) > 0):
        plotStackedHistRates(propertyLists['intraStackRatesAcceptor'], propertyLists['interStackRatesAcceptor'], ['Intra-Stack', 'Inter-Stack'], 'Acceptor', outputDir + '/18_AcceptorHoppingRate_Stacks.pdf')
        plotStackedHistTIs(propertyLists['intraStackTIsAcceptor'], propertyLists['interStackTIsAcceptor'], ['Intra-Stack', 'Inter-Stack'], 'Acceptor', outputDir + '/14_AcceptorTransferIntegral_Stacks.pdf')
    # Donor Mol Plots:
    if (len(propertyLists['intraMolRatesDonor']) > 0) or (len(propertyLists['interMolRatesDonor']) > 0):
        plotStackedHistRates(propertyLists['intraMolRatesDonor'], propertyLists['interMolRatesDonor'], ['Intra-Mol', 'Inter-Mol'], 'Donor', outputDir + '/15_DonorHoppingRate_Mols.pdf')
        plotStackedHistTIs(propertyLists['intraMolTIsDonor'], propertyLists['interMolTIsDonor'], ['Intra-Mol', 'Inter-Mol'], 'Donor', outputDir + '/11_DonorTransferIntegral_Mols.pdf')
    # Acceptor Mol Plots:
    if (len(propertyLists['intraMolRatesAcceptor']) > 0) or (len(propertyLists['interMolRatesAcceptor']) > 0):
        plotStackedHistRates(propertyLists['intraMolRatesAcceptor'], propertyLists['interMolRatesAcceptor'], ['Intra-Mol', 'Inter-Mol'], 'Acceptor', outputDir + '/17_AcceptorHoppingRate_Mols.pdf')
        plotStackedHistTIs(propertyLists['intraMolTIsAcceptor'], propertyLists['interMolTIsAcceptor'], ['Intra-Mol', 'Inter-Mol'], 'Acceptor', outputDir + '/13_AcceptorTransferIntegral_Mols.pdf')
    # Update the dataDict
    for material in chromoSpecies:
        for hopType in hopTypes:
            for hopTarget in hopTargets:
                numberOfHops = len(propertyLists[hopType + hopTarget + "Rates" + material])
                if numberOfHops == 0:
                    continue
                otherHopType = hopTypes[int((hopTypes.index(hopType) * -1) + 1)]
                proportion = numberOfHops / (numberOfHops + len(propertyLists[otherHopType + hopTarget + "Rates" + material]))
                dataDict[material.lower() + '_' + hopType + '_' + hopTarget.lower() + "hops"] = numberOfHops
                dataDict[material.lower() + '_' + hopType + '_' + hopTarget.lower() + "hops"] = proportion
    return dataDict


def plotStackedHistRates(data1, data2, labels, dataType, fileName):
    plt.figure()
    (n, bins, patches) = plt.hist([data1, data2], bins = np.logspace(1, 18, 40), stacked = True, color = ['r', 'b'], label = labels)
    plt.ylabel('Frequency (Arb. U.)')
    plt.xlabel(dataType + ' Hopping Rate (s' + r'$^{-1}$' + ')')
    plt.xlim([1,1E18])
    plt.xticks([1E0, 1E3, 1E6, 1E9, 1E12, 1E15, 1E18])
    plt.ylim([0, np.max(n) * 1.02])
    plt.legend(loc = 0, prop = {'size':18})
    plt.gca().set_xscale('log')
    plt.savefig(fileName)
    plt.close()
    print("Figure saved as", fileName)


def plotStackedHistTIs(data1, data2, labels, dataType, fileName):
    plt.figure()
    (n, bins, patches) = plt.hist([data1, data2], bins = np.linspace(0, 1.2, 20), stacked = True, color = ['r', 'b'], label = labels)
    plt.ylabel('Frequency (Arb. U.)')
    plt.xlabel(dataType + ' Transfer Integral (eV)')
    plt.xlim([0, 1.2])
    plt.ylim([0, np.max(n) * 1.02])
    plt.legend(loc = 0, prop = {'size':18})
    plt.savefig(fileName)
    plt.close()
    print("Figure saved as", fileName)


def writeCSV(dataDict, directory):
    CSVFileName = directory + '/results.csv'
    with open(CSVFileName, 'w+') as CSVFile:
        CSVWriter = csv.writer(CSVFile)
        for key, val in dataDict.items():
            CSVWriter.writerow([key, val])
    print("CSV file written to " + CSVFileName)


def createResultsPickle(directory):
    coresList = []
    for core in glob.glob(directory + '/KMC/KMClog_*.log'):
        coresList.append(re.findall(directory + '/KMC/KMClog_(.*).log', core)[0])
    keepList = []
    for core in coresList:
        selectList = []
        slot1 = directory + '/KMC/KMCslot1Results_%02d.pickle' % (int(core))
        slot2 = directory + '/KMC/KMCslot2Results_%02d.pickle' % (int(core))
        if os.path.getsize(slot1) >= os.path.getsize(slot2):
            keepList.append(slot1)
        else:
            keepList.append(slot2)
    resultsPicklesList = []
    for keeper in zip(coresList, keepList):
        newName = directory + '/KMC/KMCResults_' + str(keeper[0]) + '.pickle'
        shutil.copyfile(str(keeper[1]), newName)
        resultsPicklesList.append(newName)
    combineResultsPickles(directory, resultsPicklesList)


def combineResultsPickles(directory, pickleFiles):
    combinedData = {}
    pickleFiles = sorted(pickleFiles)
    print("%d pickle files found to combine!" % (len(pickleFiles)))
    for fileName in pickleFiles:
        # The pickle was repeatedly dumped to, in order to save time.
        # Each dump stream is self-contained, so iteratively unpickle to add the new data.
        with open(fileName, 'rb') as pickleFile:
            pickledData = pickle.load(pickleFile)
            for key, val in pickledData.items():
                try:
                    if val is None:
                        continue
                    if key not in combinedData:
                        combinedData[key] = val
                    else:
                        combinedData[key] += val
                except AttributeError:
                    pass
    # Write out the combined data
    print("Writing out the combined pickle file...")
    with open(directory + '/KMC/KMCResults.pickle', 'wb+') as pickleFile:
        pickle.dump(combinedData, pickleFile)
    print("Complete data written to", directory + "/KMCResults.pickle.")


def calculateMobility(directory, currentCarrierType, carrierData, simDims, plot3DGraphs):
    print("Considering the transport of", currentCarrierType + "...")
    print("Obtaining mean squared displacements...")
    carrierHistory, times, MSDs, timeStandardErrors, MSDStandardErrors = getCarrierData(carrierData)
    print("MSDs obtained")
    # Create the first figure that will be replotted each time
    plt.figure()
    anisotropy = plotAnisotropy(carrierData, directory, simDims, currentCarrierType, plot3DGraphs)
    if (carrierHistory is not None) and plot3DGraphs:
        print("Determining carrier hopping connections...")
        #plotConnections(chromophoreList, simDims, carrierHistory, directory, currentCarrierType)
    times, MSDs = helperFunctions.parallelSort(times, MSDs)
    print("Calculating MSD...")
    mobility, mobError, rSquared = plotMSD(times, MSDs, timeStandardErrors, MSDStandardErrors, directory, currentCarrierType)
    print("----------====================----------")
    print(currentCarrierType, "mobility for", directory, "= %.2E +- %.2E cm^{2} V^{-1} s^{-1}" % (mobility, mobError))
    print("----------====================----------")
    return mobility, mobError, rSquared, anisotropy


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--threeD", action="store_true", required=False, help="If present, use matplotlib to plot the 3D graphs (3D network, anisotropy and stack positions. This takes a while (usually a couple of minutes) to plot. Defaults to False.")
    parser.add_argument("-p", "--periodicStacks", action="store_true", required=False, help="If present, allow periodic connections to add chromophores to stacks, as well as non-periodic connections (this usually reduces the number of stacks in the system). Defaults to False.")
    parser.add_argument("-d", "--cutOffDonor", type=float, default=None, required=False, help="Specify a manual cut-off for the determination of stacks within the donor material. Connections with separation > cut-off will be classed as inter-stack. If omitted, stack cut-off will be determined automatically as the first minimum of the RDF.")
    parser.add_argument("-a", "--cutOffAcceptor", type=float, default=None, required=False, help="Specify a manual cut-off for the determination of stacks within the acceptor material. Connections with separation > cut-off will be classed as inter-stack. If omitted, stack cut-off will be determined automatically as the first minimum of the RDF.")
    parser.add_argument("-s", "--sequence", type=lambda s: [float(item) for item in s.split(',')], default=None, required=False, help='Create a figure in the current directory that describes the evolution of the anisotropy/mobility using the specified comma-delimited string as the sequence of x values. For instance -s "1.5,1.75,2.0,2.25,2.5" will assign each of the 5 following directories these x-values when plotting the mobility evolution.')
    parser.add_argument("-x", "--xlabel", default="Temperature (Arb. U.)", required=False, help='Specify an x-label for the combined plot (only used if -s is specified). Default = "Temperature (Arb. U.)"')
    args, directoryList = parser.parse_known_args()

    sys.setrecursionlimit(5000)
    holeMobilityData = []
    holeAnisotropyData = []
    electronMobilityData = []
    electronAnisotropyData = []
    dataDictList = []
    for directory in directoryList:
        # Create the figures directory if it doesn't already exist
        os.makedirs(directory + '/figures', exist_ok=True)
        # Now create the data dictionary
        dataDict = generateDataDict()
        print("\n")
        print("Getting carrier data...")
        carrierData = loadKMCResultsPickle(directory)
        print("Carrier Data obtained")
        # Now need to split up the carrierData into both electrons and holes
        carrierDataHoles, carrierDataElectrons = splitCarriersByType(carrierData)
        print("Loading chromophoreList...")
        AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = helperFunctions.loadPickle('./' + directory + '/code/' + directory + '.pickle')
        print("ChromophoreList obtained")
        morphologyShape = np.array([AAMorphologyDict[axis] for axis in ['lx', 'ly', 'lz']])
        simDims = [[-AAMorphologyDict[axis] / 2.0, AAMorphologyDict[axis] / 2.0] for axis in ['lx', 'ly', 'lz']]
        # Calculate the mobilities
        completeCarrierTypes = []
        completeCarrierData = []
        if (carrierDataHoles is not None) and (len(carrierDataHoles['ID']) > 0):
            completeCarrierTypes.append('Hole')
            completeCarrierData.append(carrierDataHoles)
        if (carrierDataElectrons is not None) and (len(carrierDataElectrons['ID']) > 0):
            completeCarrierTypes.append('Electron')
            completeCarrierData.append(carrierDataElectrons)
        for carrierTypeIndex, carrierData in enumerate(completeCarrierData):
            currentCarrierType = completeCarrierTypes[carrierTypeIndex]
            mobility, mobError, rSquared, anisotropy = calculateMobility(directory, currentCarrierType, carrierData, simDims, args.threeD)
            if currentCarrierType == 'Hole':
                holeAnisotropyData.append(anisotropy)
                holeMobilityData.append([mobility, mobError])
            elif currentCarrierType == 'Electron':
                electronAnisotropyData.append(anisotropy)
                electronMobilityData.append([mobility, mobError])
            dataDict['name'] = directory
            dataDict[currentCarrierType.lower() + '_anisotropy'] = anisotropy
            dataDict[currentCarrierType.lower() + '_mobility'] = mobility
            dataDict[currentCarrierType.lower() + '_mobility_rSquared'] = rSquared
        # Now plot the distributions!
        tempDir = directory + '/figures'
        CGToMolID = determineMoleculeIDs(CGToAAIDMaster, AAMorphologyDict, parameterDict, chromophoreList)
        dataDict = plotEnergyLevels(tempDir, chromophoreList, dataDict)
        if (args.cutOffDonor is None) or (args.cutOffAcceptor is None):
            print("No cut-off manually specified for either the donor or acceptor material, therefore automatically finding the relevant cutOff as the midpoint between the first maxmimum and the first minimum of the neighbour distance distribution...")
            print("Considering periodic neighbours is", args.periodicStacks)
            [calculatedCutOffDonor, calculatedCutOffAcceptor] = getNeighbourCutOff(chromophoreList, morphologyShape, tempDir, periodic=args.periodicStacks)
        if args.cutOffDonor is None:
            cutOffDonor = calculatedCutOffDonor
        else:
            cutOffDonor = args.cutOffDonor
        if args.cutOffAcceptor is None:
            cutOffAcceptor = calculatedCutOffAcceptor
        else:
            cutOffAcceptor = args.cutOffAcceptor
        print("Cut off in Angstroems (Donor) =", cutOffDonor)
        print("Cut off in Angstroems (Acceptor) =", cutOffAcceptor)
        stackDicts = getStacks(chromophoreList, morphologyShape, cutOffDonor, cutOffAcceptor, periodic=args.periodicStacks)
        if args.threeD:
            plotStacks3D(tempDir, chromophoreList, stackDicts, simDims)
        dataDict = plotMixedHoppingRates(tempDir, chromophoreList, parameterDict, stackDicts, CGToMolID, dataDict)
        print("\n")
        print("Writing CSV Output File...")
        writeCSV(dataDict, directory)
    print("Plotting Mobility and Anisotropy progressions...")
    if args.sequence is not None:
        if len(holeAnisotropyData) > 0:
            plotTemperatureProgression(args.sequence, holeMobilityData, holeAnisotropyData, 'Hole', args.xlabel)
        if len(electronAnisotropyData) > 0:
            plotTemperatureProgression(args.sequence, electronMobilityData, electronAnisotropyData, 'Electron', args.xlabel)
    else:
        print("Skipping plotting mobility evolution.")
