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
    os.system("mkdir -p " + morphologyDir + "/code")
    os.system("mkdir -p " + morphologyDir + "/KMC")
    os.system("mkdir -p " + morphologyDir + "/morphology")

def scpFiles(morphologyDir, kesDir):
    slashList = findIndex(morphologyDir, '/')
    morphologyName = morphologyDir[slashList[-1]+1:]
    os.system("scp kestrel:"+kesDir+"/outputFiles/"+morphologyName+"/morphology/*.log "+morphologyDir+"/morphology/")
    os.system("scp kestrel:"+kesDir+"/outputFiles/"+morphologyName+"/morphology/*.xml "+morphologyDir+"/morphology/")
    os.system("scp kestrel:"+kesDir+"/outputFiles/"+morphologyName+"/code/* "+morphologyDir+"/code/")
    os.system("scp kestrel:"+kesDir+"/outputFiles/"+morphologyName+"/KMC/KMCResults*.pickle "+morphologyDir+"/KMC")

#    os.system("scp kestrel:"+kesDir+"/outputFiles/"+morphologyName+"/KMC/KMCResults.pickle "+morphologyDir+"/KMC")

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
