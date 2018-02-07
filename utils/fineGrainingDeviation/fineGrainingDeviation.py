import os
import sys
sys.path.append('../../code')
import helperFunctions
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    AAMorphologyDir = sys.argv[1]
    AACodeDir = os.listdir(AAMorphologyDir + '/code')
    for fileName in AACodeDir:
        if ('.pickle' in fileName) and (AAMorphologyDir[AAMorphologyDir.rfind('/') + 1:] in fileName):
            AAPickleFile = AAMorphologyDir + '/code/' + fileName
            break
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = helperFunctions.loadPickle(AAPickleFile)

    boxDims = [AAMorphologyDict[x] for x in ['lx', 'ly', 'lz']]
    assert(np.linalg.norm(np.array([CGMorphologyDict[x] for x in ['lx', 'ly', 'lz']]) - np.array(boxDims)) < 0.01)

    ## DEBUG: Fix for old pickle files
    ## The old pickle AAIDs are really messed up, so this is a fix to try and get useful data from them.
    #CGChromophoreLocs = [CGMorphologyDict['position'][x] for x, atomType in enumerate(CGMorphologyDict['type']) if atomType == 'A']

    deviations = []
    for chromoID, chromophore in enumerate(chromophoreList):
        assert(chromoID == chromophore.ID)
        print("\r" + str(chromoID + 1), "of", len(chromophoreList), end = ' ')
        CGListOfPosns = []
        CGListOfMasses = []
        for CGID in chromophore.CGIDs:
            if parameterDict['CGSiteSpecies'][CGMorphologyDict['type'][CGID]] != 'None':
                CGListOfPosns.append(CGMorphologyDict['position'][CGID])
                CGListOfMasses.append(CGMorphologyDict['mass'][CGID])
        CGCoM = helperFunctions.calcCOM(CGListOfPosns, listOfMasses=CGListOfMasses)
        print(chromophore.CGIDs, chromophore.AAIDs)
        print(CGCoM, AACoM)
        input()
        deltaVec = np.array(CGCoM) - np.array(AACoM)
        for index in range(len(deltaVec)):
            if deltaVec[index] > boxDims[index] / 2.0:
                deltaVec[index] -= boxDims[index]
            elif deltaVec[index] < -boxDims[index] / 2.0:
                deltaVec[index] += boxDims[index]
        deviation = np.linalg.norm(deltaVec)
        deviations.append(deviation)
        ## DEBUG: Fix for old pickle files
        #AACoM = chromophore.posn
        #separations = []
        #for chromoLoc in CGChromophoreLocs:
        #    separation = helperFunctions.calculateSeparation(AACoM, chromoLoc)
        #    separations.append(separation)
        #deviations.append(np.min(separations))
    plt.figure()
    plt.hist(deviations, np.linspace(0,7,20), color='b')
    plt.xlim([0,7])
    plt.xlabel("CoM Deviation from CG (Ang)")
    plt.ylabel("Frequency (Arb. U.)")
    plt.savefig(AAMorphologyDir + '/figures/CoMDeviation.pdf')
    print("\nCoM Deviation saved as", AAMorphologyDir + '/figures/CoMDeviation.pdf')

    #helperFunctions.writeMorphologyXML(CGMorphologyDict, 'CG.xml')
    #helperFunctions.writeMorphologyXML(AAMorphologyDict, 'AA.xml')






    # I did this wrong: the following code checks the Chromophore Posn against the AAID CoMs, which are always identical
    #rigidBodies = {}
    #for atomID, body in enumerate(AAMorphologyDict['body']):
    #    if body != -1:
    #        if body not in rigidBodies.keys():
    #            rigidBodies[body] = [atomID]
    #        else:
    #            rigidBodies[body].append(atomID)
    ##print(rigidBodies)
    ##for rigidBody, setOfAAIDs in rigidBodies.items():
    ##    if len(set(setOfAAIDs) - set(chromophoreList[rigidBody].AAIDs)) > 0:
    ##        print("ERROR")
    ##        print(rigidBody, chromophoreList[rigidBody].ID, chromophoreList[rigidBody].AAIDs, setOfAAIDs)
    ##        input()
    #for rigidBody, setOfAAIDs in rigidBodies.items():
    #    AAIDTypes = [AAMorphologyDict['type'][AAID] for AAID in setOfAAIDs]
    #    AAIDPosns = [AAMorphologyDict['position'][AAID] for AAID in setOfAAIDs]
    #    AAIDCoM = helperFunctions.calcCOM(AAIDPosns, listOfAtomTypes=AAIDTypes)
    #    print(AAIDCoM, chromophoreList[rigidBody].posn)
    #    exit()

