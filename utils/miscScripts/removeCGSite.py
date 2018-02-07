import sys
sys.path.append('../../code')
import helperFunctions

for temperature in ['T1.0', 'T1.5', 'T1.75', 'T2.0', 'T2.25', 'T2.5']:
    fileName = './p1-L15-f0.3-P0.1-' + temperature + '-e0.1.xml'
    dictionary = helperFunctions.loadMorphologyXML(fileName, sigma=3.0)
    removeTheseAtoms = []
    for index, atomType in enumerate(dictionary['type']):
        if atomType == 'D':
            removeTheseAtoms.append(index)
    removeTheseAtoms.sort(reverse=True)
    for popID in removeTheseAtoms:
        for key in ['body', 'diameter', 'image', 'charge', 'mass', 'position', 'type']:
            dictionary[key].pop(popID)
    popTheseConstraints = [[], [], [], []]
    constraintTypes = ['bond', 'angle', 'dihedral', 'improper']
    for constraintType, key in enumerate(constraintTypes):
        for constraintID, constraint in enumerate(dictionary[key]):
            for atomID in constraint[1:]:
                if atomID in removeTheseAtoms:
                    popTheseConstraints[constraintType].append(constraintID)
                    break
    for constraintType, constraints in enumerate(popTheseConstraints):
        constraints.sort(reverse=True)
        for constraintID in constraints:
            dictionary[constraintTypes[constraintType]].pop(constraintID)
    dictionary['natoms'] -= len(removeTheseAtoms)
    helperFunctions.writeMorphologyXML(dictionary, fileName)
