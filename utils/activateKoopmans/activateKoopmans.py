import sys
sys.path.append('../../code')
import helperFunctions
from shutil import copyfile
import numpy as np

def fixTIDeltaE(chromophoreList):
    for chromo1 in chromophoreList:
        chromo1ID = chromo1.ID
        for neighbourIndex, chromo2ID in enumerate([neighbour[0] for neighbour in chromo1.neighbours]):
            chromo2 = chromophoreList[chromo2ID]
            chromo2ID = chromo2.ID
            reverseLoc = [neighbourData[0] for neighbourData in chromophoreList[chromo2ID].neighbours].index(chromo1ID)
            # TI equation is 1/2 * sqrt((HOMOSplitting)**2 - (DeltaEij)**2)
            # To activate Koopmans', need to fix this to 1/2 * HOMOSplitting
            nonKoopmansTI = chromo1.neighboursTI[neighbourIndex]
            # Check the TIs are the same for forwards and backwards
            assert(nonKoopmansTI == chromo2.neighboursTI[reverseLoc])
            # Also check these chromophores are the same species
            assert(chromo1.species == chromo2.species)
            if chromo1.species == 'Donor':
                chromo1E = chromo1.HOMO
                chromo2E = chromo2.HOMO
            else:
                chromo1E = chromo1.LUMO
                chromo2E = chromo2.LUMO
            oldDeltaEij = chromo2E - chromo1E
            # Now invert the TI equation
            koopmansTI = 0.5 * np.sqrt((nonKoopmansTI * 2)**2 + (oldDeltaEij**2))
            # Update the chromophoreList
            chromophoreList[chromo1ID].neighboursTI[neighbourIndex] = koopmansTI
            # And do the reverse
            chromophoreList[chromo2ID].neighboursTI[reverseLoc] = koopmansTI
            # Then, set the DeltaEs to 0
            chromophoreList[chromo1ID].neighboursDeltaE[neighbourIndex] = 0.0
            chromophoreList[chromo2ID].neighboursDeltaE[reverseLoc] = 0.0
    return chromophoreList


if __name__ == "__main__":
    pickleName = sys.argv[1]
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = helperFunctions.loadPickle(pickleName)
    print("Creating backup of original pickle file...")
    backedUp = pickleName + '.bak'
    try:
        with open(backedUp, 'r') as _:
            print("Backed up pickle file, " + backedUp + " already present! Running this program more than once on the same pickle will scramble the results and prevent them from being recovered. Please ensure that " + pickleName + " is a fresh pickle containing a chromophoreList where Koopmans' approximation has NOT been applied.")
            while True:
                continueFlag = input("Continue? (y/n): ")
                if continueFlag.lower() == 'y':
                    break
                elif continueFlag.lower() == 'n':
                    print("Exiting...")
                    exit()
                else:
                    print("Input not understood.")
    except FileNotFoundError:
        pass
    copyfile(pickleName, backedUp)
    print("Backup created at", backedUp + ".")
    newChromophoreList = fixTIDeltaE(chromophoreList)
    helperFunctions.writePickle((AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, newChromophoreList), pickleName)
