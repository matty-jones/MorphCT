import xml.etree.cElementTree as ET

def loadModelXML(fileName, typeDict):
    tree = ET.parse(fileName)
    root = tree.getroot()
    for model in root:
        if model.attrib['name'] == 'normalized':
            break
    constraintTypes = {'bond':{}, 'angle':{}, 'dihedral':{}}
    for constraint in model:
        if constraint.tag == 'bond':
            bondType = str(typeDict[constraint.attrib['typeA']])+'-'+str(typeDict[constraint.attrib['typeB']])
            constraintTypes['bond'][bondType] = [float(constraint.attrib['k']), float(float(constraint.attrib['l']))]
        elif constraint.tag == 'angle':
            angleType = str(typeDict[constraint.attrib['typeA']])+'-'+str(typeDict[constraint.attrib['typeB']])+'-'+str(typeDict[constraint.attrib['typeC']])
            constraintTypes['angle'][angleType] = [float(constraint.attrib['k']), float(float(constraint.attrib['theta']))]
        elif constraint.tag == 'dihedral':
            dihedralType = str(typeDict[constraint.attrib['typeA']])+'-'+str(typeDict[constraint.attrib['typeB']])+'-'+str(typeDict[constraint.attrib['typeC']])+'-'+str(typeDict[constraint.attrib['typeD']])
            constraintTypes['dihedral'][dihedralType] = [float(constraint.attrib['k1']), float(float(constraint.attrib['k2'])), float(float(constraint.attrib['k3'])), float(float(constraint.attrib['k4']))]
    return constraintTypes


def compareConstraints(modelConstraints):
    definedInModelButNotPar = []
    parBonds = [x[0] for x in par01.bondCoeffs]
    parAngles = [x[0] for x in par01.angleCoeffs]
    parDihedrals = [x[0] for x in par01.dihedralCoeffs]
    for constraintType in modelConstraints:
        if constraintType == 'bond':
            for bondType in modelConstraints[constraintType]:
                reverseBondType = '-'.join(bondType.split('-')[::-1])  # The reverse way round
                if (bondType in parBonds) or (reverseBondType in parBonds):
                    continue
                else:
                    definedInModelButNotPar.append(bondType)
        if constraintType == 'angle':
            for angleType in modelConstraints[constraintType]:
                reverseAngleType = '-'.join(angleType.split('-')[::-1])  # The reverse way round
                if (angleType in parAngles) or (reverseAngleType in parAngles):
                    continue
                else:
                    definedInModelButNotPar.append(angleType)
        if constraintType == 'dihedral':
            for dihedralType in modelConstraints[constraintType]:
                reverseDihedralType = '-'.join(dihedralType.split('-')[::-1])  # The reverse way round
                if (dihedralType in parDihedrals) or (reverseDihedralType in parDihedrals):
                    continue
                else:
                    definedInModelButNotPar.append(dihedralType)
    print definedInModelButNotPar
    exit()



if __name__ == "__main__":
    import par01
    typeDict = {'C!':'C1', 'CB':'C2', 'CW':'C3', 'CT':'C4', 'CS':'C5', 'CA':'C6', 'CP':'C1', 'HA':'H2', 'HC':'H1', 'O':'O1', 'OS':'O2', 'NS':'N1', 'S':'S1'}
    modelConstraints = loadModelXML('./model.xml', typeDict)
    compareConstraints(modelConstraints)

