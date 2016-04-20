import os
import sys
import time as T
import subprocess as sp

sys.path.append(os.getcwd()+'/code')
import fineGrainer
import helperFunctions
import runHoomd
import extractMol
import analyseMolecules

def getFilesList(direc):
    fileList = os.listdir(direc)
    morphologyFiles = []
    for fileName in fileList:
        if (fileName[-4:] == '.xml'):
            morphologyFiles.append(str(fileName))
    return morphologyFiles


def checkOutputDirectory(morphologyFile, outputDir, mode="<NO MODE>"):
    morphologyName = str(morphologyFile[:-4])
    outputDirContents = os.listdir(outputDir)
    if (morphologyName in outputDirContents):
        try:
            morphologyDirContents = os.listdir(outputDir+'/'+morphologyName+'/morphology')
        except OSError:
            os.makedirs(outputDir+'/'+morphologyName+'/morphology')
            os.makedirs(outputDir+'/'+morphologyName+'/molecules')
            return True, True
        for fileName in morphologyDirContents:
            if '.pickle' in fileName:
                #AUTO RUN
                #overwriteFlag = str(raw_input("OVERWRITE CURRENT "+mode+" DATA FOR "+morphologyName+"? (Y or N, default N): "))
                overwriteFlag = 'N'
                if ((overwriteFlag == 'y') or (overwriteFlag == 'Y')):
                    print "Cleaning directory ready for calculations..."
                    os.system('rm -rf '+str(outputDir)+'/'+morphologyName+'/*')
                    os.makedirs(str(outputDir)+'/'+morphologyName+'/morphology')
                    return True, True
                else:
                    print "Using current data..."
                    return False, True
    else:
        os.makedirs(str(outputDir)+'/'+morphologyName)        
        os.makedirs(str(outputDir)+'/'+morphologyName+'/morphology')
        os.makedirs(str(outputDir)+'/'+morphologyName+'/molecules')
    return True, True


if __name__ == '__main__':
    #### Put these into an input parameter file at some point along with everything else so that they are not defined here
    inputDir = 'inputCGMorphs'
    outputDir = 'outputFiles'
    currentDirContents = os.listdir(os.getcwd())
    if inputDir not in currentDirContents:
        os.makedirs(os.getcwd()+'/'+inputDir)
    if outputDir not in currentDirContents:
        os.makedirs(os.getcwd()+'/'+outputDir)
    inputDir = os.getcwd()+'/'+inputDir
    outputDir = os.getcwd()+'/'+outputDir
    
    morphologyFiles = getFilesList(inputDir)
    exitFlag = 0
    while exitFlag == 0:
        while True:
            print "\n---=== VALID MORPHOLOGY FILES ===---"
            for elementNo in range(len(morphologyFiles)):
                print str(elementNo)+"):", morphologyFiles[elementNo]
            print str(elementNo+1)+"): Exit program\n"
            # print "Valid files =", zip(datFiles, lammpstrjFiles)
            # AUTO RUNNING
            #runThisFile = raw_input("Please pick a file to run (integer, default = 0): ")
            try:
                runThisFile = str(sys.argv[1])
            except IndexError:
                runThisFile = '0'
            print "Selection =", runThisFile
            if len(runThisFile) == 0:
                runThisFile = 0
            else:
                try:
                    runThisFile = int(runThisFile)
                except:
                    print "Please enter an integer between 0 and", len(morphologyFiles)
                    continue
            if (runThisFile < 0) or (runThisFile > len(morphologyFiles)):
                print "Please enter an integer between 0 and", len(morphologyFiles)
                continue
            elif runThisFile == len(morphologyFiles):
                print "Exiting Program..."
                exitFlag = 1
                break
            break
        if exitFlag == 0:
            print "Checking for previous calculations of", str(morphologyFiles[runThisFile])+"..."
            runFG, runMD = checkOutputDirectory(morphologyFiles[runThisFile], outputDir, mode='MORPHOLOGY')
            if runFG == True:
                # Work out if this is a Jankowski/Marsh morphology by determining the number of hyphens in the filename
                hyphenLocs = helperFunctions.findIndex(morphologyFiles[runThisFile], '-')
                sigma = 1.
                if hyphenLocs != None:
                    if len(hyphenLocs) == 5:
                        print "This file looks like a Jankowski/Marsh morphology. The value of sigma has been set to 3 Angstroems for this morphology"
                        sigma = 3.
                t0 = T.time()
                print "Loading morphology from XML for FineGraining..."
                AAFileName, CGMorphologyDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize = fineGrainer.morphology(str(inputDir)+'/'+str(morphologyFiles[runThisFile]), sigma).analyseMorphology()
                t1 = T.time()
                elapsedTime = float(t1) - float(t0)
                if elapsedTime < 60:
                    timeunits = 'seconds.'
                elif elapsedTime < 3600:
                    elapsedTime /= 60.0
                    timeunits = 'minutes.'
                elif elapsedTime < 86400:
                    elapsedTime /= 3600.0
                    timeunits = 'hours.'
                else:
                    elapsedTime /= 86400.0
                    timeunits = 'days.'
                print "----------====================----------"
                print "FineGraining calculations completed in %.1f %s." % (float(elapsedTime), str(timeunits))
                print "----------====================----------"

            print "Running hoomd on morphology files..."
            if runMD == True:
                t2 = T.time()
                if 'AAMorphologyDict' in locals():
                    # Fine Graining calcs have been done so we already have the required data. No need to read in the pickle again.
                    morphologyFile, AAFileName, CGMorphologyDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize = runHoomd.execute(outputDir+'/'+morphologyFiles[runThisFile][:-4], AAFileName, CGMorphologyDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize) 
                else:
                    # Fine Graining calcs have not just been run, so we need to call runHoomd.py directly to unpickle the data we need
                    morphologyFile, AAFileName, CGMorphologyDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize = runHoomd.loadPickle(outputDir+'/'+morphologyFiles[runThisFile][:-4])
                t3 = T.time()
                elapsedTime = float(t3) - float(t2)
                if elapsedTime < 60:
                    timeunits = 'seconds.'
                elif elapsedTime < 3600:
                    elapsedTime /= 60.0
                    timeunits = 'minutes.'
                elif elapsedTime < 86400:
                    elapsedTime /= 3600.0
                    timeunits = 'hours.'
                else:
                    elapsedTime /= 86400.0
                    timeunits = 'days.'
                print "----------====================----------"
                print "RunHoomd calculations completed in %.1f %s." % (float(elapsedTime), str(timeunits))
                print "----------====================----------"

            print "Sorting morphology into individual molecules and outputing .xyz for DFT..."
            runMol = True
            if runMol == True:
                t4 = T.time()
                if 'AAMorphologyDict' in locals():
                    morphologyFile, AAFileName, CGMorphologyDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize = extractMol.execute(outputDir+'/'+morphologyFiles[runThisFile][:-4], AAFileName, CGMorphologyDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize)
                else:
                    morphologyFile, AAFileName, CGMorphologyDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize = extractMol.loadPickle(outputDir+'/'+morphologyFiles[runThisFile][:-4])
                t5 = T.time()
                elapsedTime = float(t5) - float(t4)
                if elapsedTime < 60:
                    timeunits = 'seconds.'
                elif elapsedTime < 3600:
                    elapsedTime /= 60.0
                    timeunits = 'minutes.'
                elif elapsedTime < 86400:
                    elapsedTime /= 3600.0
                    timeunits = 'hours.'
                else:
                    elapsedTime /= 86400.0
                    timeunits = 'days.'
                print "----------====================----------"
                print "ExtractMol calculations completed in %.1f %s." % (float(elapsedTime), str(timeunits))
                print "----------====================----------"


            print "Analysing morphology to obtain chromophores..."
       
            runAnalyse = True
            if runAnalyse == True:
                t6 = T.time()
                if 'AAMorphologyDict' in locals():
                    morphologyFile, AAFileName, CGMorphologyDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize = analyseMolecules.execute(outputDir+'/'+morphologyFiles[runThisFile][:-4], AAFileName, CGMorphologyDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize)
                else:
                    morphologyFile, AAFileName, CGMorphologyDict, AAMorphologyDict, CGtoAAIDs, moleculeAAIDs, boxSize = analyseMolecules.loadPickle(outputDir+'/'+morphologyFiles[runThisFile][:-4])
                t7 = T.time()
                elapsedTime = float(t7) - float(t6)
                if elapsedTime < 60:
                    timeunits = 'seconds.'
                elif elapsedTime < 3600:
                    elapsedTime /= 60.0
                    timeunits = 'minutes.'
                elif elapsedTime < 86400:
                    elapsedTime /= 3600.0
                    timeunits = 'hours.'
                else:
                    elapsedTime /= 86400.0
                    timeunits = 'days.'
                print "----------====================----------"
                print "AnalyseMolecules calculations completed in %.1f %s." % (float(elapsedTime), str(timeunits))
                print "----------====================----------"

            runORCA = True
            if runORCA == True:
                t8 = T.time()
                os.system('hoomd ./code/executeOrca.py '+outputDir+'/'+morphologyFiles[runThisFile][:-4])
                t9 = T.time()
                elapsedTime = float(t9) - float(t8)
                if elapsedTime < 60:
                    timeunits = 'seconds.'
                elif elapsedTime < 3600:
                    elapsedTime /= 60.0
                    timeunits = 'minutes.'
                elif elapsedTime < 86400:
                    elapsedTime /= 3600.0
                    timeunits = 'hours.'
                else:
                    elapsedTime /= 86400.0
                    timeunits = 'days.'
                print "----------====================----------"
                print "executeORCA calculations completed in %.1f %s." % (float(elapsedTime), str(timeunits))
                print "----------====================----------"

