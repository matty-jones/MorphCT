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


def updateParameters(pairData, bondCoefficients, angleCoefficients, dihedralCoefficients):
    for pairNo, pair in enumerate(ljCoeffs):
        atom1 = pair[0].split('-')[0]
        atom2 = pair[0].split('-')[1]
        epsilon = np.sqrt(pairData[atom1][0] * pairData[atom2][0])
        sigma = np.sqrt(pairData[atom1][1] * pairData[atom2][1])
        ljCoeffs[pairNo][1] = epsilon
        ljCoeffs[pairNo][2] = sigma
        dpdCoeffs[pairNo][1] = 100.0
        dpdCoeffs[pairNo][2] = 1.0
    bondCoeffs = []
    angleCoeffs = []
    dihedralCoeffs = []
    for bond in sorted(bondCoefficients.keys()):
        bondCoeffs.append([bond] + bondCoefficients[bond])
    for angle in sorted(angleCoefficients.keys()):
        angleCoeffs.append([angle] + angleCoefficients[angle])
    for dihedral in sorted(dihedralCoefficients.keys()):
        dihedralCoeffs.append([dihedral] + dihedralCoefficients[dihedral])
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
    fileName = 'par06'
    exec("from "+fileName+" import *")
    #atomRefDict = {'S1':'S', 'C1':'C!', 'C3':'CB', 'C2':'C!', 'C4':'CB', 'C5':'CW', 'C6':'CW', 'O1':'O', 'C6':'CW', 'O2':'O', 'N1':'NS', 'C8':'CT', 'C9':'CT', 'C10':'CT', 'C11':'CT', 'C12':'CT', 'C13':'CT', 'C14':'CT', 'C15':'CT', 'C16':'C!', 'C18':'CS', 'C19':'CB', 'C17':'CB', 'S2':'S', 'C20':'CA', 'C21':'CB', 'C22':'CB', 'C23':'CA', 'C24':'CS', 'C25':'CP', 'S3':'S', 'O3':'OS', 'O4':'OS' 'H1':'HC', 'H2':'HA', 'C26':'CT', 'C27':'CT', 'C28':'CT', 'C29':'CT', 'C30':'CT', 'C31':'CT', 'C32':'CT', 'C33':'CT', 'C34':'CT', 'C35':'CT', 'C36':'CT', 'C37':'CT', 'C38':'CT', 'C39':'CT', 'C40':'CT', 'C41':'CT'}
    atomRefDict = {'C!': 'C1', 'CB': 'C2', 'CW': 'C3', 'CT': 'C4', 'CS': 'C5', 'CA': 'C6', 'CP': 'C7', 'HA': 'H2', 'HC': 'H1', 'O': 'O1', 'OS': 'O2', 'NS': 'N1', 'S': 'S1'}
    pairData, bondCoefficients, angleCoefficients, dihedralCoefficients = openXMLFile('./model.xml', atomRefDict)
    ljCoeffs, dpdCoeffs, bondCoeffs, angleCoeffs, dihedralCoeffs = updateParameters(pairData, bondCoefficients, angleCoefficients, dihedralCoefficients)
    writeNewFile('./'+fileName+'.py', './'+fileName+'_modified.py', ljCoeffs, dpdCoeffs, bondCoeffs, angleCoeffs, dihedralCoeffs)
    # --=== CHECK FOR DUPLICATES ===--
#    for coeffList in ['ljCoeffs', 'dpdCoeffs', 'bondCoeffs', 'angleCoeffs', 'dihedralCoeffs']:
#        constraintNames = []
#        for constraint in eval(coeffList):
#            constraintNames.append(constraint[0])
#        print set([x for x in constraintNames if constraintNames.count(x) > 1])



