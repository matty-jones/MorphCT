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

elementaryCharge = 1.60217657E-19  # C
kB = 1.3806488E-23  # m^{2} kg s^{-2} K^{-1}
hbar = 1.05457173E-34 # m^{2} kg s^{-1}
temperature = 290  # K


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
    fileName = '3d' + carrierType + '.pdf'
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
    #plt.title('Mob = '+str(mobility)+' cm'+r'$^{2}$/Vs', y = 1.1)
    fileName = 'LinMSD' + carrierType + '.pdf'
    plt.savefig(directory + '/figures/' + fileName, bbox_inches='tight')
    plt.clf()
    print("Figure saved as", directory + "/figures/" + fileName)
    plt.semilogx(times, MSDs)
    plt.errorbar(times, MSDs, xerr = timeStandardErrors, yerr = MSDStandardErrors)
    plt.semilogx(fitX, fitY, 'r')
    plt.xlabel('Time (s)')
    plt.ylabel('MSD (m'+r'$^{2}$)')
    #plt.title('Mob = '+str(mobility)+' cm'+r'$^{2}$/Vs', y = 1.1)
    fileName = 'SemiLogMSD' + carrierType + '.pdf'
    plt.savefig(directory + '/figures/' + fileName, bbox_inches='tight')
    plt.clf()
    print("Figure saved as", directory + "/figures/" + fileName)
    plt.plot(times, MSDs)
    plt.errorbar(times, MSDs, xerr = timeStandardErrors, yerr = MSDStandardErrors)
    plt.plot(fitX, fitY, 'r')
    plt.xlabel('Time (s)')
    plt.ylabel('MSD (m'+r'$^{2}$)')
    plt.xscale('log')
    plt.yscale('log')
    #plt.title('Mob = '+str(mobility)+' cm'+r'$^{2}$/Vs', y = 1.1)
    fileName = 'LogMSD' + carrierType + '.pdf'
    plt.savefig(directory + '/figures/' + fileName, bbox_inches='tight')
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


def plotAnisotropy(carrierData, directory, simDims, carrierType):
    simExtent = [value[1] - value[0] for value in simDims]
    fig = plt.gcf()
    ax = p3.Axes3D(fig)
    xvals = []
    yvals = []
    zvals = []
    colours = []
    simDimsnm = list(map(list, np.array(simDims) / 10.))
    for carrierNo, posn in enumerate(carrierData['finalPosition']):
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
    print("----------====================----------")
    print(carrierType + " charge transport anisotropy calculated as", anisotropy)
    print("----------====================----------")
    # Reduce number of plot markers
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
    plt.savefig(directory + '/figures/anisotropy' + carrierType + '.pdf', bbox_inches='tight')
    plt.clf()
    print("Figure saved as", directory + "/figures/anisotropy" + carrierType + ".pdf")
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
    plt.savefig(fileName, bbox_inches='tight')
    plt.clf()
    print("Figure saved as " + fileName)

    plt.plot(tempData, anisotropyData, c = 'r')
    fileName = './anisotropy' + carrierType + '.pdf'
    plt.xlabel(xLabel)
    plt.ylabel(r'$\kappa$'+', Arb. U')
    plt.savefig(fileName, bbox_inches='tight')
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
        return None, None, None, None
    return binEdges, fitArgs, mean, std


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
        plt.hist(yvals, np.linspace(0,1.0,20), color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Donor Transfer Integral (eV)')
        plt.xlim([0.0, 1.2])
        fileName = 'DonorTI.pdf'

    elif mode == 'AcceptorTI':
        plt.hist(yvals, np.linspace(0,1.0,20), color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Acceptor Transfer Integral (eV)')
        plt.xlim([0.0, 1.2])
        fileName = 'AcceptorTI.pdf'

    elif mode == 'DonorTITrimmed':
        plt.hist(yvals, np.linspace(0,1.0,20), color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Donor Non-Zero Transfer Integral (eV)')
        plt.xlim([0.0, 1.2])
        fileName = 'DonorTITrimmed.pdf'

    elif mode == 'AcceptorTITrimmed':
        plt.hist(yvals, np.linspace(0,1.0,20), color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Acceptor Non-Zero Transfer Integral (eV)')
        plt.xlim([0.0, 1.2])
        fileName = 'AcceptorTITrimmed.pdf'

    elif mode == 'DonorDeltaEij':
        n, bins, patches = plt.hist(yvals, np.linspace(-0.5,0.5,20), color = ['b'])
        if gaussBins is not None:
            gaussY = gaussian(gaussBins[:-1], *fitArgs)
            scaleFactor = max(n)/max(gaussY)
            plt.plot(gaussBins[:-1], gaussY*scaleFactor, 'ro:')
        plt.ylabel('Frequency')
        plt.xlabel('Donor Delta Eij (eV)')
        plt.xlim([-0.5, 0.5])
        fileName = 'DonorDeltaEij.pdf'

    elif mode == 'AcceptorDeltaEij':
        n, bins, patches = plt.hist(yvals, np.linspace(-0.5,0.5,20), color = ['b'])
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
            plt.hist(yvals, np.linspace(0,1.0,20), color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Donor Intra-Chain TI (eV)')
        plt.xlim([0, 1.2])
        fileName = 'DonorIntraTij.pdf'

    elif mode == 'DonorInterChainTI':
        if len(yvals) > 0:
            plt.hist(yvals, np.linspace(0,1.0,20), color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Donor Inter-Chain TI (eV)')
        plt.xlim([0, 1.2])
        fileName = 'DonorInterTij.pdf'

    elif mode == 'DonorIntraChainTITrim':
        if len(yvals) > 0:
            plt.hist(yvals, np.linspace(0,1.0,20), color = ['b'])
        plt.ylabel('Frequency')
        plt.xlabel('Donor Intra-Chain TI (eV)')
        plt.xlim([0, 1.2])
        fileName = 'DonorIntraTijTrim.pdf'

    elif mode == 'DonorInterChainTITrim':
        if len(yvals) > 0:
            plt.hist(yvals, np.linspace(0,1.0,20), color = ['b'])
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

    plt.savefig(saveDir + '/' + fileName)
    plt.clf()
    print("Figure saved as", saveDir + "/" + fileName)


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



def writeCSV(dataDict, directory):
    CSVFileName = directory + '/results.csv'
    with open(CSVFileName, 'w+') as CSVFile:
        CSVWriter = csv.writer(CSVFile)
        for key, val in dataDict.items():
            CSVWriter.writerow([key, val])
    print("CSV file written to " + CSVFileName)


if __name__ == "__main__":
    sys.path.append('../../code')
    sys.path.append('../code')
    if len(sys.argv) == 1:
        for directory in os.listdir(os.getcwd()):
            if ('py' not in directory) and ('pdf' not in directory) and ('store' not in directory) and ('Store' not in directory):
                directoryList.append(directory)
    else:
        directoryList = sys.argv[1:]
    tempData = []
    holeMobilityData = []
    holeAnisotropyData = []
    electronMobilityData = []
    electronAnisotropyData = []
    combinedPlots = True
    dataDictList = []
    for directory in directoryList:
        # Create the figures directory if it doesn't already exist
        os.makedirs(directory + '/figures', exist_ok=True)
        # Now create the data dictionary
        dataDict = OrderedDict([('name', '---'), ('density', '---'), ('hole_anisotropy', '---'), ('hole_mobility', '---'), ('hole_mobility_rSquared', '---'), ('electron_anisotropy', '---'), ('electron_mobility', '---'), ('electron_mobility_rSquared', '---'),
                ('donor_HOMO_mean', '---'), ('donor_HOMO_error', '---'), ('donor_bandgap_mean', '---'), ('donor_bandgap_error', '---'), ('donor_delta_Eij_mean', '---'), ('donor_delta_Eij_std', '---'), ('donor_intra_chain_hops', '---'), ('donor_inter_chain_hops', '---'), ('donor_intra_chain_hop_proportion', '---'),
                ('acceptor_LUMO_mean', '---'), ('acceptor_LUMO_error', '---'), ('acceptor_bandgap_mean', '---'), ('acceptor_bandgap_error', '---'), ('acceptor_delta_Eij_mean', '---'), ('acceptor_delta_Eij_std', '---'), ('acceptor_intra_chain_hops', '---'), ('acceptor_inter_chain_hops', '---'), ('acceptor_intra_chain_hop_proportion', '---')])
        print("\n")
        try:
            tempData.append(getTempVal(directory))
            tempXLabel = 'T, Arb. U'
        except:
            try:
                tempData.append(getFrameVal(directory))
                tempXLabel = r'$\tau$' + ', Arb. U'
            except:
                print("No temp or frame data found in morphology name, skipping combined plots")
                combinedPlots = False
        try:
            with open(directory + '/KMC/KMCResults.pickle', 'rb') as pickleFile:
                carrierData = pickle.load(pickleFile)
        except UnicodeDecodeError:
            with open(directory + '/KMC/KMCResults.pickle', 'rb') as pickleFile:
                carrierData = pickle.load(pickleFile, encoding='latin1')
        except:
            print(sys.exc_info()[0])
            continue
        print("Carrier Data obtained")
        # Now need to split up the carrierData into both electrons and holes
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

        print("Loading chromophoreList...")
        AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = helperFunctions.loadPickle('./' + directory + '/code/' + directory + '.pickle')
        print("ChromophoreList obtained")
        simDims = [[-AAMorphologyDict[axis] / 2.0, AAMorphologyDict[axis] / 2.0] for axis in ['lx', 'ly', 'lz']]
#### NOW DO ALL OF THE BELOW BUT FOR ELECTRONS AND HOLES SEPARATELY
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
            print("Considering the transport of", currentCarrierType + "...")
            print("Obtaining mean squared displacements...")
            carrierHistory, times, MSDs, timeStandardErrors, MSDStandardErrors = getCarrierData(carrierData)
            print("MSDs obtained")
            # Create the first figure that will be replotted each time
            plt.figure()
            anisotropy = plotAnisotropy(carrierData, directory, simDims, currentCarrierType)
            if carrierHistory is not None:
                print("Determining carrier hopping connections...")
                plotConnections(chromophoreList, simDims, carrierHistory, directory, currentCarrierType)
            times, MSDs = helperFunctions.parallelSort(times, MSDs)
            print("Calculating MSD...")
            mobility, mobError, rSquared = plotMSD(times, MSDs, timeStandardErrors, MSDStandardErrors, directory, currentCarrierType)
            print("----------====================----------")
            print(currentCarrierType, "mobility for", directory, "= %.2E +- %.2E cm^{2} V^{-1} s^{-1}" % (mobility, mobError))
            print("----------====================----------")
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
        # Now we can do the plotTI/plotStacks stuff!
        # Lazy fix
        tempDir = directory + '/figures'
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

        donorBinEdges, donorFitArgs, donorMean, donorSTD = gaussFit(donorDeltaEij)
        dataDict['donor_delta_Eij_mean'] = donorMean
        dataDict['donor_delta_Eij_std'] = donorSTD
        print("Donor Delta Eij stats: mean =", donorMean, "std =", donorSTD)
        acceptorBinEdges, acceptorFitArgs, acceptorMean, acceptorSTD = gaussFit(acceptorDeltaEij)
        dataDict['acceptor_delta_Eij_mean'] = acceptorMean
        dataDict['acceptor_delta_Eij_std'] = acceptorSTD
        print("Acceptor Delta Eij stats: mean =", acceptorMean, "std =", acceptorSTD)
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
            intraChainPropn = len(donorIntraChainHop) / (len(donorIntraChainHop) + len(donorInterChainHop))
            print("Donor intra-chain hop proportion =", intraChainPropn)
            dataDict['donor_intra_chain_hops'] = len(donorIntraChainHop)
            dataDict['donor_inter_chain_hops'] = len(donorInterChainHop)
            dataDict['donor_intra_chain_hop_proportion'] = intraChainPropn
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
            intraChainPropn = len(donorIntraChainHop) / (len(donorIntraChainHop) + len(donorInterChainHop))
            print("Acceptor intra-chain hop proportion =", intraChainPropn)
            dataDict['acceptor_intra_chain_hops'] = len(acceptorIntraChainHop)
            dataDict['acceptor_inter_chain_hops'] = len(acceptorInterChainHop)
            dataDict['acceptor_intra_chain_hop_proportion'] = intraChainPropn
            plotHist(tempDir, acceptorIntraChainHop, 'AcceptorHopMix', xvals=acceptorInterChainHop)


        if len(HOMOLevels) > 0:
            print("DONOR HOMO LEVEL =", np.average(HOMOLevels), "+/-", np.std(HOMOLevels)/np.sqrt(len(HOMOLevels)))
            print("DONOR BANDGAP =", np.average(donorBandgap), "+/-", np.std(donorBandgap)/np.sqrt(len(donorBandgap)))
            dataDict['donor_HOMO_mean'] = np.average(HOMOLevels)
            dataDict['donor_HOMO_error'] = np.std(HOMOLevels)/np.sqrt(len(HOMOLevels))
            dataDict['donor_bandgap_mean'] = np.average(donorBandgap)
            dataDict['donor_bandgap_error'] = np.std(donorBandgap)/np.sqrt(len(donorBandgap))
        if len(LUMOLevels) > 0:
            print("ACCEPTOR LUMO LEVEL =", np.average(LUMOLevels), "+/-", np.std(LUMOLevels)/np.sqrt(len(LUMOLevels)))
            print("ACCEPTOR BANDGAP =", np.average(acceptorBandgap), "+/-", np.std(acceptorBandgap)/np.sqrt(len(acceptorBandgap)))
            dataDict['acceptor_HOMO_mean'] = np.average(LUMOLevels)
            dataDict['acceptor_HOMO_error'] = np.std(LUMOLevels)/np.sqrt(len(LUMOLevels))
            dataDict['acceptor_bandgap_mean'] = np.average(acceptorBandgap)
            dataDict['acceptor_bandgap_error'] = np.std(acceptorBandgap)/np.sqrt(len(acceptorBandgap))
        print("\n")
        print("Writing CSV Output File...")
        writeCSV(dataDict, directory)
    print("Plotting Mobility and Anisotropy progressions...")
    if combinedPlots is True:
        if len(holeAnisotropyData) > 0:
            plotTemperatureProgression(tempData, holeMobilityData, holeAnisotropyData, 'Hole', tempXLabel)
        if len(electronAnisotropyData) > 0:
            plotTemperatureProgression(tempData, electronMobilityData, electronAnisotropyData, 'Electron', tempXLabel)
    else:
        print("Progression plots not possible (probably due to no temperature specified). Cancelling...")
