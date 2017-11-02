import numpy as np
import xml.etree.ElementTree as ET

def openXMLFile(fileName, atomRefDict):
    modelName = 'normalized'
    pairData = {}
    bondCoefficients = {}
    angleCoefficients = {}
    dihedralCoefficients = {}
    # Load up all of the model.xml data and parse it into a useful format: ljCoeffs, bondCoeffs, angleCoeffs, dihedralCoeffs.
    # Note that we need every possible combination of ljCoeffs, but what we're looking for is a dictionary with the required data so that we can just search for par01.py['bondCoeffs']['C1-C3'] and immediately get [1.0, 1.0] to plop directly into the parameter file
    with open(fileName, 'r') as fileHandle:
        XMLFile = ET.parse(fileHandle)
    models = XMLFile.getroot()
    for model in models:
        if model.attrib['name'] != modelName:
            continue
        for constraint in model:
            if constraint.tag == 'atom':
                # LJ Pair
                newKey = atomRefDict[constraint.attrib['type']]
                pairData[newKey] = [float(constraint.attrib['epsilon']), float(constraint.attrib['diameter'])]
            elif constraint.tag == 'bond':
                # Bond
                newKey = atomRefDict[constraint.attrib['typeA']]+'-'+atomRefDict[constraint.attrib['typeB']]
                bondCoefficients[newKey] = [float(constraint.attrib['k']), float(constraint.attrib['l'])]
            elif constraint.tag == 'angle':
                # Angle
                newKey = atomRefDict[constraint.attrib['typeA']]+'-'+atomRefDict[constraint.attrib['typeB']]+'-'+atomRefDict[constraint.attrib['typeC']]
                angleCoefficients[newKey] = [float(constraint.attrib['k']), float(constraint.attrib['theta'])]
            elif constraint.tag == 'dihedral':
                # Dihedral
                newKey = atomRefDict[constraint.attrib['typeA']]+'-'+atomRefDict[constraint.attrib['typeB']]+'-'+atomRefDict[constraint.attrib['typeC']]+'-'+atomRefDict[constraint.attrib['typeD']]
                dihedralCoefficients[newKey] = [float(constraint.attrib['k1']), float(constraint.attrib['k2']), float(constraint.attrib['k3']), float(constraint.attrib['k4'])]
    # Sort the pairs by doing the geometric average
    return pairData, bondCoefficients, angleCoefficients, dihedralCoefficients


def updateParameters(sigmaScaling, epsilonScaling):
    for pairNo, pair in enumerate(ljCoeffs):
        ljCoeffs[pairNo][1] = float("%.5f" % (pair[1] * epsilonScaling))
        ljCoeffs[pairNo][2] = float("%.5f" % (pair[2] * sigmaScaling))
    for pairNo, pair in enumerate(dpdCoeffs):
        dpdCoeffs[pairNo][1] = float("%.5f" % (pair[1]))
        dpdCoeffs[pairNo][2] = float("%.5f" % (pair[2]))
    for bondNo, bond in enumerate(bondCoeffs):
        bondCoeffs[bondNo][1] = float("%.5f" % (bond[1] * epsilonScaling/(sigmaScaling**2)))
        bondCoeffs[bondNo][2] = float("%.5f" % (bond[2] * sigmaScaling))
    for angleNo, angle in enumerate(angleCoeffs):
        angleCoeffs[angleNo][1] = float("%.5f" % (angle[1] * epsilonScaling))
        angleCoeffs[angleNo][2] = float("%.5f" % (angle[2]))
    for dihedralNo, dihedral in enumerate(dihedralCoeffs):
        dihedralCoeffs[dihedralNo][1] = float("%.5f" % (dihedral[1] * epsilonScaling))
        dihedralCoeffs[dihedralNo][2] = float("%.5f" % (dihedral[2] * epsilonScaling))
        dihedralCoeffs[dihedralNo][3] = float("%.5f" % (dihedral[3] * epsilonScaling))
        dihedralCoeffs[dihedralNo][4] = float("%.5f" % (dihedral[4] * epsilonScaling))
    return ljCoeffs, dpdCoeffs, bondCoeffs, angleCoeffs, dihedralCoeffs


def writeNewFile(inputFile, outputFile, ljCoeffs, dpdCoeffs, bondCoeffs, angleCoeffs, dihedralCoeffs):
    with open(inputFile, 'r') as originalFile:
        originalLines = originalFile.readlines()
    recording = False
    lineNumbers = []
    importantParameters = ['ljCoeffs', 'dpdCoeffs', 'bondCoeffs', 'angleCoeffs', 'dihedralCoeffs']
    for lineNo, line in enumerate(originalLines):
        if recording is False:
            for importantParam in importantParameters:
                if importantParam in line:
                    lineNumbers.append([lineNo + 1])
                    recording = True
                    break
        else:
            if line == ']\n':
                lineNumbers[-1].append(lineNo)
                recording = False
    importantParameters.reverse()
    lineNumbers.reverse()
    for index, importantParam in enumerate(importantParameters):
        replacementLines = ""
        coeffList = eval(importantParam)
        for coeff in coeffList:
            replacementLines += str(coeff)+",\\\n"
        originalLines[lineNumbers[index][0]:lineNumbers[index][1]] = [replacementLines]
    with open(outputFile, 'w+') as newFile:
        newFile.writelines(originalLines)


if __name__ == "__main__":
    fileName = 'par01'
    exec("from "+fileName+" import *")
    sigmaScaling = 3.55
    epsilonScaling = 0.25
    ljCoeffs, dpdCoeffs, bondCoeffs, angleCoeffs, dihedralCoeffs = updateParameters(sigmaScaling, epsilonScaling)
    writeNewFile('./'+fileName+'.py', './'+fileName+'_modified.py', ljCoeffs, dpdCoeffs, bondCoeffs, angleCoeffs, dihedralCoeffs)
    # --=== CHECK FOR DUPLICATES ===--
#    for coeffList in ['ljCoeffs', 'dpdCoeffs', 'bondCoeffs', 'angleCoeffs', 'dihedralCoeffs']:
#        constraintNames = []
#        for constraint in eval(coeffList):
#            constraintNames.append(constraint[0])
#        print set([x for x in constraintNames if constraintNames.count(x) > 1])



