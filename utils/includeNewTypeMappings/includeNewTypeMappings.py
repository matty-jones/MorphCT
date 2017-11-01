import sys
sys.path.append('../../code')
import helperFunctions
import fineGrainer

if __name__ == "__main__":
    pickleFile = sys.argv[1]
    AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList = helperFunctions.loadPickle(pickleFile)
    morphologyXML = parameterDict['inputMorphDir'] + '/' + parameterDict['morphology']
    system = fineGrainer.morphology(morphologyXML, parameterDict['morphology'][:-4], parameterDict, chromophoreList)
    print("Obtaining the newTypeMappings...")
    newTypeMappings = system.getNewTypeMappings(system.CGToTemplateDirs, system.CGToTemplateForceFields)
    print("Updating parameter dictionary with the newTypeMappings...")
    parameterDict['newTypeMappings'] = newTypeMappings
    helperFunctions.writePickle((AAMorphologyDict, CGMorphologyDict, CGToAAIDMaster, parameterDict, chromophoreList), pickleFile)
