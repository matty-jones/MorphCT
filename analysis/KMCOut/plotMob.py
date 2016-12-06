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


def combinePickleFiles(directory):
    pickleFiles = []
    for fileName in os.listdir(directory):
        if 'KMCData' == fileName[:7]:
            pickleFiles.append(directory + '/' + fileName)
    print len(pickleFiles), "pickle files found to concatenate."
    print "Concatenating data..."
    carrierList = []
    for pickleNo, pickleFileName in enumerate(pickleFiles):
        print "\rLoading data", pickleNo + 1, "of", len(pickleFiles),
        sys.stdout.flush()
        with open(pickleFileName, 'r') as pickleFile:
            carrierList += pickle.load(pickleFile)
    print "\nAll data concatenated!"
    print "Writing combined pickle..."
    with open(directory + '/combinedKMCData.pickle', 'w+') as pickleFile:
        pickle.dump(carrierList, pickleFile)
    return carrierList


def getData(carrierList):
    squaredDisps = {}
    actualTimes = {}
    noChromophoresVisited = {}
    completeCarrierHistory = lil_matrix(carrierList[0].carrierHistoryMatrix.shape, dtype = int)
    totalDataPoints = 0
    totalDataPointsAveragedOver = 0
    for carrier in carrierList:
        if (carrier.currentTime > carrier.lifetime * 2) or (carrier.currentTime < carrier.lifetime / 2.0) or (carrier.noHops == 1):
            totalDataPoints += 1
            continue
        carrierKey = str(carrier.lifetime)
        if carrierKey not in squaredDisps:
            squaredDisps[carrierKey] = [(carrier.displacement * 1E-10) ** 2]  # Carrier displacement is in angstroems, convert to metres
            actualTimes[carrierKey] = [carrier.currentTime]
            noChromophoresVisited[carrierKey] = [len(findNonZero(carrier.carrierHistoryMatrix)[0])]
        else:
            squaredDisps[carrierKey].append((carrier.displacement * 1E-10) ** 2)  # Carrier displacement is in angstroems, convert to metres
            actualTimes[carrierKey].append(carrier.currentTime)
            noChromophoresVisited[carrierKey].append(len(findNonZero(carrier.carrierHistoryMatrix)[0]))
        completeCarrierHistory += carrier.carrierHistoryMatrix
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
    return completeCarrierHistory, times, MSDs, timeStandardErrors, MSDStandardErrors


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
    plt.savefig('./3d.pdf')
    plt.show()
    exit()
    pass


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
    fileName = 'LinMSD.png'
    plt.savefig(directory + '/' + fileName)
    plt.clf()
    print "Figure saved as", directory + "/" + fileName
    plt.semilogx(times, MSDs)
    plt.errorbar(times, MSDs, xerr = timeStandardErrors, yerr = MSDStandardErrors)
    plt.semilogx(fitX, fitY, 'r')
    plt.xlabel('Time (s)')
    plt.ylabel('MSD (m'+r'$^{2}$)')
    #plt.title('Mob = '+str(mobility)+' cm'+r'$^{2}$/Vs', y = 1.1)
    fileName = 'LogMSD.png'
    plt.savefig(directory + '/' + fileName)
    plt.clf()
    print "Figure saved as", directory + "/" + fileName
    return mobility, mobError


if __name__ == "__main__":
    sys.path.append('../../code')
    directory = os.getcwd() + '/' + sys.argv[1]
    print "Combining Pickle Files..."
    needToCombine = True
    for fileName in os.listdir(directory):
        if 'combinedKMCData' in fileName:
            combinedPickleName = fileName
            needToCombine = False
            break
    if needToCombine is True:
        print "combinedKMCData.pickle not found! Combining data..."
        carrierList = combinePickleFiles(directory)
    else:
        with open(directory + '/combinedKMCData.pickle', 'r') as pickleFile:
            carrierList = pickle.load(pickleFile)
    print "Carrier List obtained"
    print "Obtaining mean squared displacements..."
    carrierHistory, times, MSDs, timeStandardErrors, MSDStandardErrors = getData(carrierList)
    print "MSDs obtained"
    # Create the first figure that will be replotted each time
    plt.figure()
    #plotHeatMap(carrierHistory, directory)
    # READ IN THE MAIN CHROMOPHORELIST PICKLE FILE TO DO THIS
    print "Loading chromophoreList..."
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList = helperFunctions.loadPickle('./' + sys.argv[1] + '/' + sys.argv[1] + '.pickle')
    print "ChromophoreList obtained"
    plotConnections(chromophoreList, [AAMorphologyDict['lx'], AAMorphologyDict['ly'], AAMorphologyDict['lz']], carrierHistory, directory)
    times, MSDs = helperFunctions.parallelSort(times, MSDs)
    mobility, mobError = plotMSD(times, MSDs, timeStandardErrors, MSDStandardErrors, directory)
    print "----------====================----------"
    print "Mobility for", directory[helperFunctions.findIndex(directory, '/')[-1] + 1:], "= %.2E +- %.2E cm^{2} V^{-1} s^{-1}" % (mobility, mobError)
    print "----------====================----------"
