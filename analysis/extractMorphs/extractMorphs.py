import sys
sys.path.append('../../../code')
import helperFunctions

if __name__ == "__main__":
    temperature = sys.argv[1]
    fileName = 'p1-L15-f0.0-P0.1-T' + temperature + '-e0.5'
    pickleFile = '../../' + fileName + '/code/' + fileName + '.pickle'
    data = helperFunctions.loadPickle(pickleFile)
    helperFunctions.writeMorphologyXML(data[0], './' + fileName + '_AA.xml')
