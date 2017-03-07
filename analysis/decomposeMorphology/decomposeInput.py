import sys
sys.path.append('../../code')
import helperFunctions


if __name__ == "__main__":
    morphFile = sys.argv[1]
    CGMorphology = helperFunctions.loadMorphologyXML('./' + morphFile)
    fullereneIndices = []
    siteTypes = ['C']
    for index, siteType in enumerate(CGMorphology['type']):
        if siteType in siteTypes:
            fullereneIndices.append(index)
    popList = sorted(fullereneIndices, reverse=True)
    print "There are", len(popList), "atoms to remove that are type", str(siteTypes) + "."
    print "Removing atoms..."
    properties = ['position', 'image', 'mass', 'diameter', 'type', 'body', 'charge']
    for atomProperty in properties:
        for index in popList:
            CGMorphology[atomProperty].pop(index)
    CGMorphology['natoms'] -= len(popList)
    print "Done! Updating constraints..."
    constraints = ['bond', 'angle', 'dihedral', 'improper']
    for constraintName in constraints:
        constraintsToRemove = []
        for constraintID, constraint in enumerate(CGMorphology[constraintName]):
            for atomIndex in constraint[1:]:
                if atomIndex in popList:
                    constraintsToRemove.append(constraintID)
                    break
        print "Found", len(constraintsToRemove), constraintName, "constraints to remove..."
        popList = sorted(constraintsToRemove, reverse=True)
        for index in popList:
            CGMorphology[constraintName].pop(index)
    print "Done! Writing XML..."
    helperFunctions.writeMorphologyXML(CGMorphology, './TRIM' + morphFile)
