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


def getNNCutOff(chromophoreList, outputDir):
    separationDist = []
    for chromo1 in chromophoreList:
        for chromo2 in chromophoreList:
            if chromo1.ID == chromo2.ID:
                continue
            separation = np.linalg.norm(chromo2.posn - chromo1.posn)
            if separation < 10:
                separationDist.append(separation)
    plt.figure()
    (n, binEdges, patches) = plt.hist(separationDist, bins = 20)
    plt.savefig(outputDir + "/neighbourHist.pdf")
    plt.close()
    print("Neighbour histogram figure saved as", outputDir + "/neighbourHist.pdf")
    bins = 0.5*(binEdges[1:]+binEdges[:-1])
    bins = np.insert(bins, 0, 0)
    n = np.insert(n, 0, 0)
    nonZero = False
    minimumFound = False
    count = 0
    lastNonZeroIndex = 0
    nextNonZeroIndex = 0
    for index, number in enumerate(n):
        if count == 2:
            minimumFound = True
        if nonZero is True:
            if number == 0:
                count += 1
            else:
                if minimumFound is True:
                    nextNonZeroIndex = index
                    break
                lastNonZeroIndex = index
                count = 0
        if number != 0:
            nonZero = True
    print(bins, n)
    print("Minimum between", bins[lastNonZeroIndex], bins[nextNonZeroIndex])
    cutOff = 0.5 * (bins[lastNonZeroIndex] + bins[nextNonZeroIndex])
    return cutOff


def getStacks(chromophoreList, morphologyShape, cutOff):
    # Create a neighbourlist based on the cutoff
    neighbourDict = createNeighbourList(chromophoreList, morphologyShape, cutOff)
    # Do the usual stackList neighbourList stuff
    stackList = [_ for _ in range(len(chromophoreList))]
    for stackID in range(len(stackList)):
        stackList = updateStack(stackID, stackList, neighbourDict)
    print("There are", len(set(stackList)), "stacks in the system")
    stackDict = {}
    for index, chromophore in enumerate(chromophoreList):
        stackDict[chromophore.ID] = stackList[index]
    return stackDict


def createNeighbourList(chromophoreList, morphologyShape, cutOff):
    neighbourDict = {}
    for chromo1 in chromophoreList:
        for [chromo2ID, relImage] in chromo1.neighbours:
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


def plot3DMorphology(outputDir, chromophoreList, stackDict):
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
    ax.set_xlim([-30, 30])
    ax.set_ylim([-30, 30])
    ax.set_zlim([-30, 30])
    plt.savefig(outputDir + "/stacks.pdf")
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
            rate = calculateHopRate(lambdaij * elementaryCharge, Tij * elementaryCharge, deltaE * elementaryCharge, T)
            if stackDict[chromo.ID] == stackDict[chromo.neighbours[index][0]]:
                intraStackRates.append(rate)
            else:
                interStackRates.append(rate)
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
    tempDirs = []
    for fileName in os.listdir(os.getcwd()):
        if ('py' not in fileName) and ('pdf' not in fileName) and ('store' not in fileName):
            tempDirs.append(fileName)
    for tempDir in tempDirs:
        print("\n")
        for fileName in os.listdir(os.getcwd() + '/' + tempDir):
            if ("pickle" in fileName) and (tempDir in fileName):
                mainMorphologyPickleName = os.getcwd() + '/' + tempDir + '/' + fileName
        AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = helperFunctions.loadPickle(mainMorphologyPickleName)
        # Check chromo coherence first
        for index, chromophore in enumerate(chromophoreList):
            assert index == chromophore.ID, "index != chromophore.ID: %r != %r" % (index, chromophore.ID)
        cutOff = getNNCutOff(chromophoreList, tempDir)
        print("Cut off in Angstroems =", cutOff)
        morphologyShape = np.array([AAMorphologyDict[_] for _ in ['lx', 'ly', 'lz']])
        stackDict = getStacks(chromophoreList, morphologyShape, cutOff)
        plot3DMorphology(tempDir, chromophoreList, stackDict)
        plotMixedHoppingRate(tempDir, chromophoreList, stackDict)
