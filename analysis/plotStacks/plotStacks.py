import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import csv
sys.path.append('../../code')
import helperFunctions
try:
    import mpl_toolkits.mplot3d as p3
except ImportError:
    print("Could not import 3D plotting engine, calling the plot3DMorphology function will result in an error!")


elementaryCharge = 1.60217657E-19 # C
kB = 1.3806488E-23 # m^{2} kg s^{-2} K^{-1}
hbar = 1.05457173E-34 # m^{2} kg s^{-1}


def getNeighbourCutOff(chromophoreList, morphologyShape, outputDir, periodic=True):
    separationDist = []
    for chromo1 in chromophoreList:
        for chromo2Details in chromo1.neighbours:
            if (chromo2Details is None) or ((periodic is False) and (not np.array_equal(chromo2Details[1], [0, 0, 0]))) or (chromo1.ID == chromophoreList[chromo2Details[0]].ID):
                continue
            chromo2 = chromophoreList[chromo2Details[0]]
            separation = np.linalg.norm((np.array(chromo2.posn) + (np.array(chromo2Details[1]) * np.array(morphologyShape))) - chromo1.posn)
            separationDist.append(separation)
    plt.figure()
    (n, binEdges, patches) = plt.hist(separationDist, bins = 20)
    plt.savefig(outputDir + "/neighbourHist.pdf")
    plt.close()
    print("Neighbour histogram figure saved as", outputDir + "/neighbourHist.pdf")
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
    cutOff = (bins[maximaIndices[0]] + bins[minimaIndices[0]]) / 2.0
    return cutOff


def getStacks(chromophoreList, morphologyShape, cutOff, periodic=True):
    # Create a neighbourlist based on the cutoff
    neighbourDict = createNeighbourList(chromophoreList, morphologyShape, cutOff, periodic)
    # Do the usual stackList neighbourList stuff
    stackList = [_ for _ in range(len(chromophoreList))]
    for stackID in range(len(stackList)):
        stackList = updateStack(stackID, stackList, neighbourDict)
    print("There are", len(set(stackList)), "stacks in the system")
    stackDict = {}
    for index, chromophore in enumerate(chromophoreList):
        stackDict[chromophore.ID] = stackList[index]
    return stackDict


def createNeighbourList(chromophoreList, morphologyShape, cutOff, periodic=True):
    neighbourDict = {}
    for chromo1 in chromophoreList:
        for [chromo2ID, relImage] in chromo1.neighbours:
            if periodic is False:
                if not np.array_equal(relImage, [0, 0, 0]):
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


def plot3DMorphology(outputDir, chromophoreList, stackDict, simDims):
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    colours = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'w']
    stackList = {}
    for chromophore in chromophoreList:
        stackID = stackDict[chromophore.ID]
        if stackID not in stackList.keys():
            stackList[stackID] = [chromophore]
        else:
            stackList[stackID].append(chromophore)
    for stackID, chromos in enumerate(stackList.values()):
        for chromo in chromos:
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
    plt.savefig(outputDir + "/stacks.pdf", bbox_inches='tight')
    plt.close()
    print("3D Stack figure saved as", outputDir + "/stacks.pdf")
    #plt.show()


def plotMixedHoppingRate(outputDir, chromophoreList, stackDict):
    intraStackRates = []
    interStackRates = []
    for chromo in chromophoreList:
        for index, Tij in enumerate(chromo.neighboursTI):
            deltaE = chromo.neighboursDeltaE[index]
            lambdaij = 0.130
            T = 290
            try:
                rate = calculateHopRate(lambdaij * elementaryCharge, Tij * elementaryCharge, deltaE * elementaryCharge, T)
                if stackDict[chromo.ID] == stackDict[chromo.neighbours[index][0]]:
                    intraStackRates.append(rate)
                else:
                    interStackRates.append(rate)
            except TypeError:
                pass
    print(len(intraStackRates), len(interStackRates))
    plt.figure()
    plt.hist([intraStackRates, interStackRates], bins = np.logspace(1, 18, 40), stacked = True, color = ['r', 'b'], label = ['Intra-Stack', 'Inter-Stack'])
    plt.ylabel('Frequency')
    plt.xlabel('Donor Hopping Rate (s' + r'$^{-1}$' + ')')
    plt.xlim([1,1E18])
    plt.xticks([1E0, 1E3, 1E6, 1E9, 1E12, 1E15, 1E18])
    #plt.ylim([0,8000])
    plt.legend(loc = 2, prop = {'size':18})
    plt.gca().set_xscale('log')
    fileName = outputDir + '/DonorHoppingRateMixed.pdf'
    plt.savefig(fileName)
    plt.close()


def calculateHopRate(lambdaij, Tij, deltaEij, T):
    # Semiclassical Marcus Hopping Rate Equation
    kij = ((2 * np.pi) / hbar) * (Tij ** 2) * np.sqrt(1.0 / (4 * lambdaij * np.pi * kB * T)) * np.exp(-((deltaEij + lambdaij)**2) / (4 * lambdaij * kB * T))
    return kij


if __name__ == "__main__":
    periodic = True
    try:
        cutOff = float(sys.argv[1])
        tempDirs = sys.argv[2:]
    except ValueError:
        cutOff = None
        tempDirs = sys.argv[1:]
    sys.setrecursionlimit(5000)
    for tempDir in tempDirs:
        print("\n")
        for fileName in os.listdir(os.getcwd() + '/' + tempDir):
            if ("pickle" in fileName) and (tempDir in fileName):
                mainMorphologyPickleName = os.getcwd() + '/' + tempDir + '/' + fileName
        AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = helperFunctions.loadPickle(mainMorphologyPickleName)
        morphologyShape = np.array([AAMorphologyDict[axis] for axis in ['lx', 'ly', 'lz']])
        simDims = [[-axis/2.0, axis/2.0] for axis in morphologyShape]
        # Check chromo coherence first
        for index, chromophore in enumerate(chromophoreList):
            assert index == chromophore.ID, "index != chromophore.ID: %r != %r" % (index, chromophore.ID)
        if cutOff is None:
            print("No cut-off manually specified, therefore automatically finding cutOff as the midpoint between the first maxmimum and the first minimum of the neighbour distance distribution.")
            print("Considering periodic neighbours is", periodic)
            cutOff = getNeighbourCutOff(chromophoreList, morphologyShape, tempDir, periodic=periodic)
        print("Cut off in Angstroems =", cutOff)
        stackDict = getStacks(chromophoreList, morphologyShape, cutOff, periodic=periodic)
        plot3DMorphology(tempDir, chromophoreList, stackDict, simDims)
        plotMixedHoppingRate(tempDir, chromophoreList, stackDict)
