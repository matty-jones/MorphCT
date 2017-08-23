import sys
sys.path.append('../../code')

import helperFunctions

if __name__ == "__main__":
    fileName = sys.argv[1]
    template = helperFunctions.loadMorphologyXML(fileName)
    for constraintName in ['angle','dihedral','improper']:
        for constraint in template[constraintName]:
            for index in range(1, len(constraint) - 1):
                shouldBeBond = [constraint[index], constraint[index + 1]]
                correct = False
                for bond in template['bond']:
                    if (bond[1:] == shouldBeBond) or (bond[-1:0:-1] == shouldBeBond): # Bond in reverse
                        correct = True
                        break
                if correct is False:
                    raise SystemError("PROBLEM. The " + constraintName + " " + repr(constraint) + " should be between the bonded atoms " + repr(shouldBeBond) + ", however this bond does not exist in the template.")
    print "All constraints occur between bonded atoms, so no problem there."
