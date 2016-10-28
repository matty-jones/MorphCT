class FFData:
    def __init__(self, dataFile):
        with open(dataFile, 'r') as dataFileHandle:
            dataLines = dataFileHandle.readlines()
        record = False
        print len(dataLines[0])
        for lineNo, line in enumerate(dataLines):
            if len(line) == 0:
                record = False
                continue
            props = ['atom', 'vdw', 'bond', 'angle', 'torsion', 'charge']
            for propertyName in props:
                if propertyName in line:
                    print filter(None, line.split(' '))
                    exit()





if __name__ == "__main__":
    OPLSAA = FFData('./oplsaa.prm.txt')
