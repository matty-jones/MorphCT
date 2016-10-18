import os
import sys


def getMorphologyDirs(directory):
    morphologyDirs = []
    parentDir = os.listdir(directory)
    for dirName in parentDir:
        if "-T" in dirName:
            morphologyDirs.append(directory+"/"+dirName)
    return morphologyDirs


def makeDirectoryStructure(morphologyDir, kesDir):
    currentStructure = os.listdir(morphologyDir)
    if "M15" not in kesDir:
        if "M01T0" not in currentStructure:
            os.system("mkdir "+morphologyDir+"/M01TI0")
            os.system("mkdir "+morphologyDir+"/M01TI0/KMC")
        else:
            inDir = os.listdir(morphologyDir+"/M01TI0")
            if "KMC" not in inDir:
                os.system("mkdir "+morphologyDir+"/M01TI0/KMC")
    else:
        if "M15T0" not in currentStructure:
            os.system("mkdir "+morphologyDir+"/M15TI0")
            os.system("mkdir "+morphologyDir+"/M15TI0/KMC")
        else:
            inDir = os.listdir(morphologyDir+"/M15TI0")
            if "KMC" not in inDir:
                os.system("mkdir "+morphologyDir+"/M15TI0/KMC")


def scpFiles(morphologyDir, kesDir):
    slashList = findIndex(morphologyDir, '/')
    morphologyName = morphologyDir[slashList[-1]+1:]
    if "M15" not in kesDir:
        os.system("scp kestrel:"+kesDir+"/outputFiles/"+morphologyName+"/morphology/*.pickle "+morphologyDir+"/")
        os.system("scp kestrel:"+kesDir+"/outputFiles/"+morphologyName+"/morphology/*.xml "+morphologyDir+"/")
        os.system("scp kestrel:"+kesDir+"/outputFiles/"+morphologyName+"/chromophores/*.csv "+morphologyDir+"/M01TI0")
        os.system("scp kestrel:"+kesDir+"/outputFiles/"+morphologyName+"/KMC/*.csv "+morphologyDir+"/M01TI0/KMC")
    else:
        os.system("scp kestrel:"+kesDir+"/outputFiles/"+morphologyName+"/morphology/*.pickle "+morphologyDir+"/")
        os.system("scp kestrel:"+kesDir+"/outputFiles/"+morphologyName+"/morphology/*.xml "+morphologyDir+"/")
        os.system("scp kestrel:"+kesDir+"/outputFiles/"+morphologyName+"/chromophores/*.csv "+morphologyDir+"/M15TI0")
        os.system("scp kestrel:"+kesDir+"/outputFiles/"+morphologyName+"/KMC/*.csv "+morphologyDir+"/M15TI0/KMC")


def findIndex(string, character):
    '''This function returns the locations of an inputted character in an inputted string'''
    index = 0
    locations = []
    while index < len(string):
        if string[index] == character:
            locations.append(index)
        index += 1
    if len(locations) == 0:
        return None
    return locations


if __name__ == "__main__":
    kesDir = sys.argv[1]
    directory = os.getcwd()
    morphologyDirs = getMorphologyDirs(directory)
    for morphologyDir in morphologyDirs:
        makeDirectoryStructure(morphologyDir, kesDir)
        scpFiles(morphologyDir, kesDir)
