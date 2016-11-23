import os
import sys
import numpy as np
import random as R
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import helperFunctions
try:
    import mpl_toolkits.mplot3d.axes3d as p3
except ImportError:
    print "Could not import 3D plotting engine, calling the plotCarrier subroutine will result in an error!"


plottingSubroutines = False
elementaryCharge = 1.60217657E-19 # C
kB = 1.3806488E-23 # m^{2} kg s^{-2} K^{-1}
hbar = 1.05457173E-34 # m^{2} kg s^{-1}

# Flag for outputting all hop times to a CSV (slow!)
debug = False
###


class carrier:
    def __init__(self, chromophoreList, parameterDict, chromoID, lifetime):
        if parameterDict['recordCarrierHistory'] is True:
            self.carrierHistoryMatrix = np.zeros((len(chromophoreList), len(chromophoreList)), dtype = int)
        else:
            self.carrierHistoryMatrix = None
        self.image = [0, 0, 0]
        self.currentChromophore = chromophoreList[chromoID]
        self.lifetime = lifetime
        self.currentTime = 0.0


def execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList):
    # Determine the maximum simulation times based on the parameter dictionary
    simulationTimes = np.logspace(np.log10(parameterDict['minimumSimulationTime']), np.log10(parameterDict['maximumSimulationTime']), len(parameterDict['procIDs']))
    carrierList = []
    # For each specifed carrier
    for carrierNo in range(len(numberOfCarriersPerSimulationTime)):
        # For each specified simulation time (do the loops this way round to even out the job lists when we split all the KMC jobs over the CPUs)
        for lifetime in simulationTimes:
            # Find an initial position (for now, just pick a chromophore at random)
            startChromoID = R.randint(0, len(chromophoreList) - 1)
            carrierList.append(carrier(chromophoreList, parameterDict, chromoID, lifetime))
    # The carrierList is now like the ORCAJobsList, so split it over each procID

    # Then run a bunch of singleCoreRunKMC jobs

    # In that, probably make a csv file with the important data, to save us having to dump to the pickle every 5 minutes


    return AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList


if __name__ == "__main__":
    try:
        pickleFile = sys.argv[1]
    except:
        print "Please specify the pickle file to load to continue the pipeline from this point."
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList = helperFunctions.loadPickle(pickleFile)
    execute(AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList, carrierList)







































































class chargeCarrier:
    def __init__(self, initialChromophore, singlesData, TIDict, boxSize, temperature, CTOutputDir):
        self.TIDict = TIDict
        self.singlesData = singlesData
        # SinglesData Dictionary is of the form:
        # key = realChromoID, val = [xPos, yPos, zPos, HOMO-1, HOMO, LUMO, LUMO+1]
        self.initialPosition = np.array(singlesData[initialChromophore][0:3])
        self.currentChromophore = initialChromophore
        self.imagePosition = [0, 0, 0]
        self.globalTime = 0
        self.T = temperature
        self.boxSize = boxSize
        self.simDims = [[-boxSize[0]/2.0, boxSize[0]/2.0], [-boxSize[1]/2.0, boxSize[1]/2.0], [-boxSize[2]/2.0, boxSize[2]/2.0]]
        self.CTOutputDir = CTOutputDir
        self.hopHistory = [int(initialChromophore)]
        self.reinitialise()


    def reinitialise(self):
        self.position = np.array(self.singlesData[self.currentChromophore][0:3])
        self.chromoLength = int(self.singlesData[self.currentChromophore][7])
        if plottingSubroutines == True:
            plotCarrier(self.singlesData, self.TIDict, self.currentChromophore, self.position, self.imagePosition, self.initialPosition, self.simDims, self.globalTime)


    def calculateLambdaij(self):
        # The equation for the internal reorganisation energy was obtained from the data given in
        # Johansson, E and Larsson, S; 2004, Synthetic Metals 144: 183-191.
        # External reorganisation energy obtained from 
        # Liu, T and Cheung, D. L. and Troisi, A; 2011, Phys. Chem. Chem. Phys. 13: 21461-21470
        lambdaExternal = 0.11 # eV
        if self.chromoLength < 12:
            lambdaInternal = 0.20826 - (self.chromoLength*0.01196)
        else:
            lambdaInternal = 0.06474
        lambdaeV = lambdaExternal+lambdaInternal
        lambdaJ = lambdaeV*elementaryCharge
        return lambdaJ


    def calculateEij(self, destination):
        Ei = self.singlesData[self.currentChromophore][4] # HOMO LEVEL
        Ej = self.singlesData[destination][4]
        deltaEijeV = Ej - Ei
        deltaEijJ = deltaEijeV*elementaryCharge
        return deltaEijJ


    def calculateHopRate(self, lambdaij, Tij, deltaEij):
        # Error in Lan 2008, should be just 1/hbar not squared
        kij = ((2*np.pi)/hbar)*(Tij**2)*np.sqrt(1.0/(4*lambdaij*np.pi*kB*self.T))*np.exp(-((deltaEij+lambdaij)**2)/(4*lambdaij*kB*self.T))
        #print "Prefactor =", ((2*np.pi)/hbar)*(Tij**2)*np.sqrt(1.0/(4*lambdaij*np.pi*kB*self.T))
        #print "Exponent =", (np.exp(-((deltaEij+lambdaij)**2)/(4*lambdaij*kB*self.T)))
        # Durham code had a different prefactor == Tij**2/hbar * sqrt(pi/(lambda*kB*T))
        return kij


    def determineHopTime(self, rate):
        if rate != 0:
            while True:
                x = R.random()
                if (x != 0.0) and (x != 1.0):
                    break
            tau = -np.log(x)/rate
        else:
            # Zero rate, therefore set the hop time to very long
            tau = 1E20
        return tau


    def calculateHop(self):
        hopTimes = []
        lambdaij = self.calculateLambdaij()
        for hopTarget in self.TIDict[self.currentChromophore]:
            #print "\n"
            transferIntegral = hopTarget[1]*elementaryCharge
            deltaEij = self.calculateEij(hopTarget[0])
            hopRate = self.calculateHopRate(lambdaij, transferIntegral, deltaEij)
            hopTime = self.determineHopTime(hopRate)
            #print "For this hop target:", hopTarget
            #print "TI =", transferIntegral
            #print "deltaEij =", deltaEij
            #print "hopRate =", hopRate
            #print "hopTime =", hopTime
            hopTimes.append([hopTarget[0], hopTime])
        hopTimes.sort(key = lambda x:x[1]) # Sort by ascending hop time
        deltaEij = self.calculateEij(hopTimes[0][0])
        #print "\n"
        #print hopTimes
        # if deltaEij <= 0.0:
        #     print "Downstream hop", deltaEij
        # else:
        #     print "Upstream hop", deltaEij
        # print "HopTime =", hopTimes[0][1]
        # raw_input('Post Hop Pause')
        self.performHop(hopTimes[0][0], hopTimes[0][1])


    def performHop(self, destinationChromophore, hopTime):
        initialPosition = np.array(self.singlesData[self.currentChromophore][0:3])
        destinationPosition = np.array(self.singlesData[destinationChromophore][0:3])
        deltaPosition = destinationPosition - initialPosition
        # print "Hopping from", self.currentChromophore, "to", destinationChromophore
        # print "Current =", initialPosition, "Destination =", destinationPosition
        # print "deltaPosition =", deltaPosition
        # Work out if we've crossed a periodic boundary
        for axis in range(3):
            if np.absolute(deltaPosition[axis]) > self.boxSize[axis]/2.0:
                # Crossed a periodic boundary condition. Find out which one!
                if destinationPosition[axis] > initialPosition[axis]:
                    # Crossed a negative boundary, so increment image
                    self.imagePosition[axis] -= 1
                elif destinationPosition[axis] < initialPosition[axis]:
                    # Crossed a positive boundary
                    self.imagePosition[axis] += 1
        # Image sorted, now move the charge
        self.currentChromophore = destinationChromophore
        # Increment the time
        self.globalTime += hopTime
        self.hopHistory.append(int(destinationChromophore))
        if debug == True:
            fileName = self.CTOutputDir+'/hopTimes.csv'
            writeDebugCSV(fileName, hopTime)
        self.reinitialise()


def plotCarrier(singleChromos, TIDict, currentChromo, currentPosn, image, initialPosition, simDims, globalTime):
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    # Draw the simulation box first
    ax.plot([simDims[0][0], simDims[0][1]], [simDims[1][0], simDims[1][0]], [simDims[2][0], simDims[2][0]], c = 'k')
    ax.plot([simDims[0][0], simDims[0][1]], [simDims[1][1], simDims[1][1]], [simDims[2][0], simDims[2][0]], c = 'k')
    ax.plot([simDims[0][0], simDims[0][1]], [simDims[1][0], simDims[1][0]], [simDims[2][1], simDims[2][1]], c = 'k')
    ax.plot([simDims[0][0], simDims[0][1]], [simDims[1][1], simDims[1][1]], [simDims[2][1], simDims[2][1]], c = 'k')

    ax.plot([simDims[0][0], simDims[0][0]], [simDims[1][0], simDims[1][1]], [simDims[2][0], simDims[2][0]], c = 'k')
    ax.plot([simDims[0][0], simDims[0][0]], [simDims[1][0], simDims[1][1]], [simDims[2][1], simDims[2][1]], c = 'k')
    ax.plot([simDims[0][1], simDims[0][1]], [simDims[1][0], simDims[1][1]], [simDims[2][0], simDims[2][0]], c = 'k')
    ax.plot([simDims[0][1], simDims[0][1]], [simDims[1][0], simDims[1][1]], [simDims[2][1], simDims[2][1]], c = 'k')

    ax.plot([simDims[0][0], simDims[0][0]], [simDims[1][0], simDims[1][0]], [simDims[2][0], simDims[2][1]], c = 'k')
    ax.plot([simDims[0][0], simDims[0][0]], [simDims[1][1], simDims[1][1]], [simDims[2][0], simDims[2][1]], c = 'k')
    ax.plot([simDims[0][1], simDims[0][1]], [simDims[1][0], simDims[1][0]], [simDims[2][0], simDims[2][1]], c = 'k')
    ax.plot([simDims[0][1], simDims[0][1]], [simDims[1][1], simDims[1][1]], [simDims[2][0], simDims[2][1]], c = 'k')

    # Draw current chromo position
    ax.scatter(currentPosn[0], currentPosn[1], currentPosn[2], s = 50, c = 'r')
    # Draw neighbouring chromos
    for hopOption in TIDict[currentChromo]:
        hopOptionPosn = np.array(singleChromos[hopOption[0]][0:3])
        ax.scatter(hopOptionPosn[0], hopOptionPosn[1], hopOptionPosn[2], s = 20, c = 'b')

    # Complete plot
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    extraSpace = 5
    ax.set_xlim(simDims[0][0]-extraSpace, simDims[0][1]+extraSpace)
    ax.set_ylim(simDims[1][0]-extraSpace, simDims[1][1]+extraSpace)
    ax.set_zlim(simDims[2][0]-extraSpace, simDims[2][1]+extraSpace)

    displacement = helperFunctions.calculateSeparation(initialPosition, currentPosn)
    displacement = list(str(displacement))
    displacementString = ''
    for i in range(5):
        try:
            displacementString += displacement[i]
        except IndexError:
            break

    plt.title(str(image)+", Disp ="+str(displacementString)+", t = "+str(globalTime))
    fileList = os.listdir('./')
    fileNameCounter = 0
    for files in fileList:
        if ".png" in files:
            fileNameCounter += 1
    fileNameAddon = str(fileNameCounter)
    while len(fileNameAddon) < 3:
        fileNameAddon = '0'+fileNameAddon
    plt.savefig('./test'+fileNameAddon+'.png')
    # plt.show()
    # raw_input('Break for Ctrl-C')
    print "Image saved as ./test"+fileNameAddon+".png"
    plt.close(fig)


def randomPosition(boxSize, singleChromos):
    randX = R.uniform(-boxSize[0]/2.0, boxSize[0]/2.0)
    randY = R.uniform(-boxSize[1]/2.0, boxSize[1]/2.0) 
    randZ = R.uniform(-boxSize[2]/2.0, boxSize[2]/2.0)
    randomPosn = [randX, randY, randZ]

    separationToSegments = [] # Of the form [ [seg1No, seg1Sep], [seg2No, seg2Sep] ... ]
    for segNo, chromo in singleChromos.iteritems():
        separationToSegments.append([segNo, helperFunctions.calculateSeparation(randomPosn, np.array(chromo[1:4]))])
    separationToSegments.sort(key = lambda x: x[1])
    return int(separationToSegments[0][0])


def writeDebugCSV(fileName, data):
    with open(fileName, 'a+') as fileHandle:
        document = csv.writer(fileHandle, delimiter = ',')
        document.writerow([data])


def writeCSVFile(fileName, carrierNo, displacement, hops, time, hopHistory):
    with open(fileName, 'a+') as fileHandle:
        document = csv.writer(fileHandle, delimiter = ',')
        displacementInM = float(displacement)*1E-10 # Convert from angstroms to metres
        document.writerow([carrierNo, displacementInM, hops, time]+hopHistory)
    print "CSV file,", str(fileName)+", appended to successfully."
    

def loadChromophoreFiles(CSVDir, CTOutputDir):
    failed = False
    singlesData = {}
    pairsData = []
    try:
        with open(CSVDir+'/singles.csv', 'r') as singlesFile:
            singlesReader = csv.reader(singlesFile, delimiter=',')
            for row in singlesReader:
                try:
                    singlesData[int(float(row[0]))] = map(float, row[1:])
                except IndexError:
                    print "Nonstandard line in singles.csv!"
                    continue
        with open(CSVDir+'/pairs.csv', 'r') as pairsFile:
            pairsReader = csv.reader(pairsFile, delimiter=',')
            for row in pairsReader:
                try:
                    pairsData.append([int(float(row[0])), int(float(row[1]))] + [x for x in map(float, row[2:])])
                except IndexError:
                    print "Nonstandard line in pairs.csv!"
                    continue
    except IOError:
        print "CSV files singles.csv and pairs.csv not found in the chromophores directory."
        print "Please run transferIntegrals.py to generate these files from the ORCA outputs."
        failed = True
    TIDict = {}
    totalPairs = 0
    numberOfZeroes = 0
    for pair in pairsData:
        if pair[0] not in TIDict:
            TIDict[pair[0]] = []
        if pair[1] not in TIDict:
            TIDict[pair[1]] = []
        TIDict[pair[0]].append([pair[1], pair[-1]]) # Neighbour and corresponding Tij
        TIDict[pair[1]].append([pair[0], pair[-1]]) # Reverse Hop
        if pair[-1] == 0.0:
            numberOfZeroes += 1
        totalPairs += 1
    print "There are", totalPairs, "total possible hop destinations, and", numberOfZeroes, "of them have a transfer integral of zero ("+str(int(float(numberOfZeroes)/float(totalPairs)*100))+"%)..."
    return singlesData, pairsData, TIDict, failed


def getPreviousCarriers(csvFileName):
    try:
        with open(csvFileName, 'r') as csvFile:
            csvLines = csv.reader(csvFile, delimiter=',')
            for row in csvLines:
                pass
        nextCarrierNo = int(row[0])+1
    except:
        nextCarrierNo = 0
    return nextCarrierNo


def runCarrier(singlesData, TIDict, boxSize, temperature, CTOutputDir):
    # Loop Start Here
    # Pick a random chromophore to inject to
    while True:
        initialChromophore = randomPosition(boxSize, singlesData)
        #print "Injecting onto", initialChromophore, "TIDict =", TIDict[initialChromophore]
        # Initialise a carrier
        hole = chargeCarrier(initialChromophore, singlesData, TIDict, boxSize, temperature, CTOutputDir)
        numberOfHops = 0
        newGlobalTime = 0.0
        # Start hopping!
        while True:
            if limitByHops == True:
                if numberOfHops == hopLimit:
                    break
                print "Performing hop number", numberOfHops+1
            else:
                if newGlobalTime > simulationTime:
                    break
            hole.calculateHop()
            newGlobalTime = hole.globalTime
            numberOfHops += 1
        if numberOfHops != 1:
            # Put a catch in in case we start on a trap site where hopping time > simulation Time
            # (otherwise this messes up the average time calculation)
            break
    initialPos = hole.initialPosition
    currentPos = np.array([hole.position[0]+(hole.imagePosition[0]*boxSize[0]), hole.position[1]+(hole.imagePosition[1]*boxSize[1]), hole.position[2]+(hole.imagePosition[2]*boxSize[2])])
    displacement = helperFunctions.calculateSeparation(initialPos, currentPos)
    hopHistory = list(set(hole.hopHistory))
    return displacement, numberOfHops, newGlobalTime, hopHistory



def execute(morphologyName, boxSize, simulationTime):
    #R.seed(32)
    CSVDir = os.getcwd()+'/outputFiles/'+morphologyName+'/chromophores'
    CTOutputDir = os.getcwd()+'/outputFiles/'+morphologyName+'/KMC'
    singlesData, pairsData, TIDict, failed = loadChromophoreFiles(CSVDir, CTOutputDir)
    if failed == True:
        return
    if limitByHops == True:
        csvFileName = CTOutputDir+'/CTHops_'+str(hopLimit)+'.csv'
    else:
        csvFileName = CTOutputDir+'/CTTime_'+str(simulationTime)+'.csv'
    # Load the CSV file for this configuration (if it exists)
    nextCarrierNo = getPreviousCarriers(csvFileName)
    firstRun = True
    if nextCarrierNo < numberOfCarriersToSimulate:
        for carrierNo in range(numberOfCarriersToSimulate):
            if firstRun == True:
                if carrierNo != nextCarrierNo:
                    continue
                else:
                    firstRun = False
            displacement, numberOfHops, newGlobalTime, hopHistory = runCarrier(singlesData, TIDict, boxSize, temperature, CTOutputDir)
            # Update CSV file
            if plottingSubroutines == False:
                print numberOfHops, "hops complete. Displacement of", displacement, "for carrier number", carrierNo, "in", newGlobalTime, "s."
                writeCSVFile(csvFileName, carrierNo, displacement, numberOfHops, newGlobalTime, hopHistory)
            else:
                print numberOfHops, "hops complete. Graphs plotted. Simulation terminating. No CSV data will be saved while plotting == True."
                break

    
    
def loadPickle(morphologyFile, simulationTime):
    morphologyName = morphologyFile[helperFunctions.findIndex(morphologyFile,'/')[-1]+1:]
    outputDir = './outputFiles'
    morphologyList = os.listdir(outputDir)
    for allMorphologies in morphologyList:
        if morphologyName in allMorphologies:
            outputDir += '/'+morphologyName
            break
    pickleFound = False
    for fileName in os.listdir(outputDir+'/morphology'):
        if fileName == morphologyName+'.pickle':
            pickleLoc = outputDir+'/morphology/'+fileName
            pickleFound = True
    if pickleFound == False:
        print "Pickle file not found. Please run morphCT.py again to create the required HOOMD inputs."
        exit()
    print "Pickle found at", str(pickleLoc)+"."
    print "Loading data..."
    with open(pickleLoc, 'r') as pickleFile:
        (AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize) = pickle.load(pickleFile)
    execute(morphologyName, boxSize, simulationTime)
    return morphologyFile, AAfileName, CGMoleculeDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize


if __name__ == "__main__":
    morphologyFile = sys.argv[1]
    try:
        simulationTime = float(sys.argv[2])
    except:
        simulationTime = float(defaultSimTime)
    loadPickle(morphologyFile, simulationTime)
