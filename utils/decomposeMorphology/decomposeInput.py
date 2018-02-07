import sys
sys.path.append('../../code')
import helperFunctions


if __name__ == "__main__":
    morphFile = sys.argv[1]
    CGMorphology = helperFunctions.loadMorphologyXML('./' + morphFile)
    fullereneIndices = []
    #allIndices = list(range(len(CGMorphology['position'])))
    #indicesToKeep = [675, 676, 677, 678, 679, 680, 681, 682, 683, 684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719, 1125, 1126, 1127, 1128, 1129, 1130, 1131, 1132, 1133, 1134, 1135, 1136, 1137, 1138, 1139, 1140, 1141, 1142, 1143, 1144, 1145, 1146, 1147, 1148, 1149, 1150, 1151, 1152, 1153, 1154, 1155, 1156, 1157, 1158, 1159, 1160, 1161, 1162, 1163, 1164, 1165, 1166, 1167, 1168, 1169]
    #fullereneIndices = list(set(allIndices) - set(indicesToKeep))
    siteTypes = ['C2', 'C3', 'C4', 'O1', 'O2', 'H1', 'H2']
    for index, siteType in enumerate(CGMorphology['type']):
        if siteType in siteTypes:
            fullereneIndices.append(index)
    popList = sorted(fullereneIndices, reverse=True)
    print("There are", len(popList), "atoms to remove that are type", str(siteTypes) + ".")
    print("Removing atoms...")
    properties = ['position', 'image', 'mass', 'diameter', 'type', 'body', 'charge']
    for atomProperty in properties:
        for index in popList:
            CGMorphology[atomProperty].pop(index)
    CGMorphology['natoms'] -= len(popList)
    print("Done! Updating constraints...")
    constraints = ['bond', 'angle', 'dihedral', 'improper']
    for constraintName in constraints:
        constraintsToRemove = []
        for constraintID, constraint in enumerate(CGMorphology[constraintName]):
            for atomIndex in constraint[1:]:
                if atomIndex in popList:
                    constraintsToRemove.append(constraintID)
                    break
        print("Found", len(constraintsToRemove), constraintName, "constraints to remove...")
        popList = sorted(constraintsToRemove, reverse=True)
        for index in popList:
            CGMorphology[constraintName].pop(index)
    print("Done! Writing XML...")
    helperFunctions.writeMorphologyXML(CGMorphology, './TRIM' + morphFile)
