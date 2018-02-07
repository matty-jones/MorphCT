import sys
sys.path.append('../../code')
import helperFunctions


if __name__ == "__main__":
    fileName = sys.argv[1]
    AAMorph = helperFunctions.loadMorphologyXML(fileName)
    rigidBodies = sorted(list(set(AAMorph['body'])))
    rigidBodies.remove(-1)
    rigidBodyLookup = {-1: -1}
    for index, bodyNo in enumerate(rigidBodies):
        rigidBodyLookup[bodyNo] = index
    for atomID, rigidBodyID in enumerate(AAMorph['body']):
        AAMorph['body'][atomID] = rigidBodyLookup[rigidBodyID]
    helperFunctions.writeMorphologyXML(AAMorph, 'bodyFix_' + fileName)
