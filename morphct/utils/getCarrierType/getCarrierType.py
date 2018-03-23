import os
import sys
import pickle
import numpy as np
from scipy.sparse import lil_matrix
sys.path.append('../../code/')
import helperFunctions


def addInCarrierTypes(carrierData, chromophoreList):
    print("Predicting carrier types given their initial position and the chromophore locations...")
    if 'carrierType' not in carrierData.keys():
        carrierData['carrierType'] = []
    for carrierIndex, initialPosition in enumerate(carrierData['initialPosition']):
        carrierType = chromophoreTypeLookupFromPosn(initialPosition, chromophoreList)
        carrierData['carrierType'].append(carrierType)
    print("Splitting the carrierHistoryMatrix into its hole and electron components...")
    # Now need to split the carrierHistoryMatrix up into a hole and an electron matrix
    matrixShape = carrierData['carrierHistoryMatrix'].shape
    holeHistoryMatrix = lil_matrix(matrixShape, dtype = int)
    electronHistoryMatrix = lil_matrix(matrixShape, dtype = int)
    nonZero = carrierData['carrierHistoryMatrix'].nonzero()
    for coordIndex, chromo1ID in enumerate(nonZero[0]):
        chromo2ID = nonZero[1][coordIndex]
        type1 = chromophoreTypeLookupFromID(chromo1ID, chromophoreList)
        type2 = chromophoreTypeLookupFromID(chromo2ID, chromophoreList)
        if type1 != type2:
            raise SystemError("Found a hop between a donor and an acceptor! Something went seriously wrong!")
        if type1 == 'Donor':
            holeHistoryMatrix[chromo1ID,chromo2ID] = carrierData['carrierHistoryMatrix'][chromo1ID,chromo2ID]
        elif type1 == 'Acceptor':
            electronHistoryMatrix[chromo1ID,chromo2ID] = carrierData['carrierHistoryMatrix'][chromo1ID,chromo2ID]
    carrierData['holeHistoryMatrix'] = holeHistoryMatrix
    carrierData['electronHistoryMatrix'] = electronHistoryMatrix
    return carrierData


def chromophoreTypeLookupFromPosn(position, chromophoreList):
    '''A function that looks up the species of a chromophore based on a position'''
    for chromophore in chromophoreList:
        if (np.array(chromophore.posn) == np.array(position)).all():
            if chromophore.species == 'Donor':
                return 'Hole'
            elif chromophore.species == 'Acceptor':
                return 'Electron'
            break
    raise SystemError("No chromophore found, or chromophore species not present! Check inputs!")


def chromophoreTypeLookupFromID(ID, chromophoreList):
    '''A function that looks up the species of a chromophore based on its ID'''
    chromophore = chromophoreList[ID]
    if chromophore.ID != ID:
        raise SystemError("Chromophore List Indices scrambled!")
    return chromophore.species


if __name__ == "__main__":
    sys.path.append('../../code')
    directoryList = []
    for directory in os.listdir(os.getcwd()):
        if ('py' not in directory) and ('pdf' not in directory) and ('store' not in directory):
            directoryList.append(directory)
    for directory in directoryList:
        print("Loading the carrier data...")
        try:
            with open(directory + '/KMCResults.pickle', 'rb') as pickleFile:
                carrierData = pickle.load(pickleFile)
        except:
            raise SystemError("No carrier data found for "+directory)
        print("Carrier Data obtained")
        print("Loading the chromophore data...")
        try:
            pickleFile = directory + '/' + directory + '.pickle'
            AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = helperFunctions.loadPickle(pickleFile)
        except:
            raise SystemError("No chromophore data found for "+directory)
        print("ChromophoreList obtained")
        carrierData = addInCarrierTypes(carrierData, chromophoreList)
        print("Updating carrier data...")
        with open(directory + '/KMCResults_SplitCarrier.pickle', 'wb+') as pickleFile:
            pickle.dump(carrierData, pickleFile)
        print("New carrier file written to", directory + '/KMCResults_SplitCarrier.pickle')
