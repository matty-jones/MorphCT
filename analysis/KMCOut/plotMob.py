import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import cPickle as pickle
import scipy.optimize
import scipy.stats
from scipy.sparse import lil_matrix
from scipy.sparse import find as findNonZero
sys.path.append('../../code/')
import helperFunctions
try:
    import mpl_toolkits.mplot3d as p3
except ImportError:
    print "Could not import 3D plotting engine, calling the plotMolecule3D function will result in an error!"

elementaryCharge = 1.60217657E-19  # C
kB = 1.3806488E-23  # m^{2} kg s^{-2} K^{-1}
temperature = 290  # K


def getData(carrierData):
    carrierHistory = carrierData['carrierHistoryMatrix']
    totalDataPoints = 0
    totalDataPointsAveragedOver = 0
    squaredDisps = {}
    actualTimes = {}
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
        totalDataPointsAveragedOver += 1
        totalDataPoints += 1
    times = []
    MSDs = []
    timeStandardErrors = []
    MSDStandardErrors = []
    for time, disps in squaredDisps.iteritems():
        times.append(float(time))
        timeStandardErrors.append(np.std(actualTimes[time]) / len(actualTimes[time]))
        MSDs.append(np.average(disps))
        MSDStandardErrors.append(np.std(disps) / len(disps))
    return carrierHistory, times, MSDs, timeStandardErrors, MSDStandardErrors


def plotHeatMap(carrierHistory, directory):
    # A simple function that displays the carrier history matrix as a heatmap
    # First convert carrierHistory from a sparse to a dense, binary matrix
    binaryHistory = np.zeros(carrierHistory.shape)
    nonZero = np.nonzero(carrierHistory)
    for i, xval in enumerate(nonZero[0]):
        binaryHistory[xval, nonZero[1][i]] = 1.0
    #carrierHistory = carrierHistory.toarray().astype(float)
    #carrierHistory /= np.amax(carrierHistory)
    fig = plt.gcf()

    ax = fig.add_subplot(111)
    plt.imshow(binaryHistory)
    ax.set_aspect('equal')

    cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    cax.patch.set_alpha(0)
    cax.set_frame_on(False)
    plt.colorbar(orientation='vertical')
    plt.savefig(directory + '/heatmap.pdf')
    plt.clf()


def plotConnections(chromophoreList, simExtent, carrierHistory, directory):
    # A complicated function that shows connections between carriers in 3D that carriers prefer to hop between.
    # Connections that are frequently used are highlighted in black, whereas rarely used connections are more white.
    fig = plt.gcf()
    ax = p3.Axes3D(fig)
    # Find a good normalisation factor
    carrierHistory = carrierHistory.toarray()
    normalizeTo = np.log(np.max(carrierHistory))
    for chromo1, row in enumerate(carrierHistory):
        for chromo2, value in enumerate(row):
            if value > 0:
                coords1 = chromophoreList[chromo1].posn
                coords2 = chromophoreList[chromo2].posn
                #ax.scatter(coords1[0], coords1[1], coords1[2], c = 'k', s = '5')
                #ax.scatter(coords2[0], coords2[1], coords2[2], c = 'k', s = '5')
                line = [coords2[0] - coords1[0], coords2[1] - coords1[1], coords2[2] - coords2[1]]
                if (np.abs(coords2[0] - coords1[0]) < simExtent[0] / 2.0) and (np.abs(coords2[1] - coords1[1]) < simExtent[1] / 2.0) and (np.abs(coords2[2] - coords1[2]) < simExtent[2] / 2.0):
                    colourIntensity = np.log(value) / normalizeTo
                    ax.plot([coords1[0], coords2[0]], [coords1[1], coords2[1]], [coords1[2], coords2[2]], c = plt.cm.jet(colourIntensity), linewidth = 0.5)
    fileName = '3d.pdf'
    plt.savefig(directory + '/' + fileName)
    print "Figure saved as", directory + "/" + fileName
    plt.clf()
    #plt.show()


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


def plotMSD(times, MSDs, timeStandardErrors, MSDStandardErrors, directory):
    fit = np.polyfit(times, MSDs, 1)
    fitX = np.linspace(np.min(times), np.max(times), 100)
    gradient, intercept, rVal, pVal, stdErr = scipy.stats.linregress(times, MSDs)
    print "StandardError", stdErr
    print "Fitting rVal =", rVal
    fitY = (fitX * gradient) + intercept
    mobility, mobError = calcMobility(fitX, fitY, np.average(timeStandardErrors), np.average(MSDStandardErrors))
    plt.plot(times, MSDs)
    plt.errorbar(times, MSDs, xerr = timeStandardErrors, yerr = MSDStandardErrors)
    plt.plot(fitX, fitY, 'r')
    plt.xlabel('Time (s)')
    plt.ylabel('MSD (m'+r'$^{2}$)')
    #plt.title('Mob = '+str(mobility)+' cm'+r'$^{2}$/Vs', y = 1.1)
    fileName = 'LinMSD.pdf'
    plt.savefig(directory + '/' + fileName)
    plt.clf()
    print "Figure saved as", directory + "/" + fileName
    plt.semilogx(times, MSDs)
    plt.errorbar(times, MSDs, xerr = timeStandardErrors, yerr = MSDStandardErrors)
    plt.semilogx(fitX, fitY, 'r')
    plt.xlabel('Time (s)')
    plt.ylabel('MSD (m'+r'$^{2}$)')
    #plt.title('Mob = '+str(mobility)+' cm'+r'$^{2}$/Vs', y = 1.1)
    fileName = 'LogMSD.pdf'
    plt.savefig(directory + '/' + fileName)
    plt.clf()
    print "Figure saved as", directory + "/" + fileName
    return mobility, mobError


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


def plotAnisotropy(carrierData, directory):
    fig = plt.gcf()
    ax = p3.Axes3D(fig)
    xvals = []
    yvals = []
    zvals = []
    for image in carrierData['image']:
        xvals.append(image[0])
        yvals.append(image[1])
        zvals.append(image[2])
    anisotropy = calculateAnisotropy(xvals, yvals, zvals)
    print "Anisotropy calculated as", anisotropy
    plt.scatter(xvals, yvals, zs = zvals, c = 'b', s = 20)
    plt.scatter(0, 0, zs = 0, c = 'r', s = 50)
    plt.savefig(directory + '/anisotropy.pdf')
    plt.clf()
    print "Figure saved as", directory + "/anisotropy.pdf"
    return anisotropy


def getTempVal(string):
    hyphenList = helperFunctions.findIndex(string, '-')
    tempVal = float(string[hyphenList[-2] + 2 : hyphenList[-1]])
    return tempVal


def plotTemperatureProgression(tempData, mobilityData, anisotropyData):
    plt.gcf()
    xvals = tempData
    yvals = list(np.array(mobilityData)[:,0])
    yerrs = list(np.array(mobilityData)[:,1])
    plt.xlabel('Temperature, Arb. U')
    plt.ylabel('Mobility, cm'+r'$^{2}$ '+'V'+r'$^{-1}$'+r's$^{-1}$')
    plt.title('p1-L15-f0.0-P0.1-TX.X-e0.1', fontsize = 24)
    plt.xlim([1.4, 2.6])
    plt.semilogy(xvals, yvals, c = 'b')
    plt.errorbar(xvals, yvals, xerr = 0, yerr = yerrs)
    fileName = './mobility.pdf'
    plt.savefig(fileName)
    plt.clf()
    print "Figure saved as " + fileName

    plt.plot(tempData, anisotropyData, c = 'r')
    fileName = './anisotropy.pdf'
    plt.xlabel('Temperature, Arb. U')
    plt.ylabel('Anisotropy, Arb. U')
    plt.savefig(fileName)
    plt.clf()
    print "Figure saved as " + fileName


if __name__ == "__main__":
    sys.path.append('../../code')
    directoryList = []
    for directory in os.listdir(os.getcwd()):
        if '-T' in directory:
            directoryList.append(directory)
    tempData = []
    mobilityData = []
    anisotropyData = []
    for directory in directoryList:
        tempData.append(getTempVal(directory))
        with open(directory + '/KMCResults.pickle', 'r') as pickleFile:
            carrierData = pickle.load(pickleFile)
        print "Carrier Data obtained"
        print "Obtaining mean squared displacements..."
        carrierHistory, times, MSDs, timeStandardErrors, MSDStandardErrors = getData(carrierData)
        print "MSDs obtained"
        # Create the first figure that will be replotted each time
        plt.figure()
        anisotropy = plotAnisotropy(carrierData, directory)
        #plotHeatMap(carrierHistory, directory)
        # READ IN THE MAIN CHROMOPHORELIST PICKLE FILE TO DO THIS
        print "Loading chromophoreList..."
        AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, emptyCarrierList = helperFunctions.loadPickle('./' + directory + '/' + directory + '.pickle')
        print "ChromophoreList obtained"
        print "Determining carrier hopping connections..."
        plotConnections(chromophoreList, [AAMorphologyDict['lx'], AAMorphologyDict['ly'], AAMorphologyDict['lz']], carrierHistory, directory)
        times, MSDs = helperFunctions.parallelSort(times, MSDs)
        print "Calculating MSD..."
        mobility, mobError = plotMSD(times, MSDs, timeStandardErrors, MSDStandardErrors, directory)
        print "----------====================----------"
        print "Mobility for", directory, "= %.2E +- %.2E cm^{2} V^{-1} s^{-1}" % (mobility, mobError)
        print "----------====================----------"
        anisotropyData.append(anisotropy)
        mobilityData.append([mobility, mobError])
    print "Plotting temperature progression..."
    plotTemperatureProgression(tempData, mobilityData, anisotropyData)
