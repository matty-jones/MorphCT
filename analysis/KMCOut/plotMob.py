import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import cPickle as pickle
import scipy.optimize
import scipy.stats
from scipy.sparse import lil_matrix
from scipy.sparse import find as findNonZero


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
    noChromophoresVisited = {}
    completeCarrierHistory = lil_matrix(carrierList[0].shape, dtype = int)
    totalDataPoints = 0
    totalDataPointsAveragedOver = 0
    for carrier in carrierList:
        if (carrier.currentTime > carrier.lifetime * 2) or (carrier.currentTime < carrier.lifetime / 2.0) or (carrier.noHops == 1):
            totalDataPoints += 1
            continue
        carrierKey = str(carrier.lifetime)
        if carrierKey not in squaredDisps:
            squaredDisps[carrierKey] = [carrier.displacement ** 2]
            noChromophoresVisited[carrierKey] = [len(findNonZero(carrier.carrierHistoryMatrix)[0])]
        else:
            squaredDisps[carrierKey].append(carrier.displacement ** 2)
            noChromophoresVisited[carrierKey].append(len(findNonZero(carrier.carrierHistoryMatrix)[0]))
        completeCarrierHistory += carrier.carrierHistoryMatrix
        totalDataPointsAveragedOver += 1
        totalDataPoints += 1
    times = []
    MSDs = []
    for time, disps in squaredDisps.iteritems():
        times.append(float(time))
        MSDs.append(np.average(disps))




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
    print "Carrier List obtained!"
    getData(carrierList)
