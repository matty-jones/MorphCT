import cPickle as pickle
import numpy as np
import sys
from scipy.sparse import lil_matrix
sys.path.append('../../../code')
import helperFunctions

def chromoIDToMol(CGToAAIDMaster, chromophoreList):
    CGIDToMol = {}
    for molID, molecule in enumerate(CGToAAIDMaster):
        CGIDs = molecule.keys()
        for CGID in CGIDs:
            CGIDToMol[CGID] = molID
    chromoIDToMolDict = {}
    for chromo in chromophoreList:
        chromoIDToMolDict[chromo.ID] = CGIDToMol[chromo.CGIDs[0]]
    return chromoIDToMolDict


if __name__ == "__main__":
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList = helperFunctions.loadPickle('p3httest.pickle')
    print "Determining the molIDs for each chromophore"
    chromoIDToMolDict = chromoIDToMol(CGToAAIDMaster, chromophoreList)
    for i in range(10):
        hopTypes = []
        hopQuantities = []
        KMCResultsPickle = 'KMCResults_%02d.pickle' % (i)
        with open('./'+KMCResultsPickle, 'r') as pickleFile:
            KMCResults = pickle.load(pickleFile)
        nonZeroData = KMCResults['carrierHistoryMatrix'].nonzero()
        for index, chromo1 in enumerate(nonZeroData[0]):
            chromo2 = nonZeroData[1][index]
            if chromo1 < chromo2:
                numberOfHops = KMCResults['carrierHistoryMatrix'][chromo1, chromo2]
                if chromoIDToMolDict[chromophoreList[chromo1].ID] == chromoIDToMolDict[chromophoreList[chromo2].ID]:
                    hopTypes.append(1)
                    for j in range(numberOfHops):
                        hopQuantities.append(1)
                else:
                    hopTypes.append(0)
                    for j in range(numberOfHops):
                        hopQuantities.append(0)
        print "For", KMCResultsPickle
        print "The number of unique intra-chain routes was", np.sum(hopTypes)
        print "The number of unique inter-chain routes was", len(hopTypes) - np.sum(hopTypes)
        print "The ratio of intra- to inter-chain routes was", np.sum(hopTypes) / float(len(hopTypes))
        print "The total number of intra-chain hops was", np.sum(hopQuantities)
        print "The total number of inter-chain hops was", len(hopQuantities) - np.sum(hopQuantities)
        print "The ratio of intra- to inter-chain hops was", np.sum(hopQuantities) / float(len(hopQuantities))
        raw_input("PAUSE...")
