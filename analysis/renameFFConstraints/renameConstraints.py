import sys
sys.path.append('../../code')
import helperFunctions

if __name__ == "__main__":
    fileName = sys.argv[1]
    constraintTypes = ['bond', 'angle', 'dihedral']
    # Create the dictionary of the template file using loadMorphologyXML
    templateDict = helperFunctions.loadMorphologyXML(fileName)
    for constraintType in constraintTypes:  # Iterate over all constraints
        for index, constraint in enumerate(templateDict[constraintType]):
            newConstraint0 = ""  # Create a new constraint label
            for atomID in constraint[1:]:
                newConstraint0 += templateDict['type'][atomID]
                newConstraint0 += '-'
                # Assign the constraint label
            templateDict[constraintType][index][0] = newConstraint0[:-1]
    # Write the morphology xml
    helperFunctions.writeMorphologyXML(templateDict, fileName)
