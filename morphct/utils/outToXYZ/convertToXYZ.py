import sys
import os

if __name__ == "__main__":
    files = os.listdir('./')
    for fileName in files:
        if '.out' not in fileName:
            continue
        print(fileName)
        with open('./'+fileName, 'r') as outFile:
            outData = outFile.readlines()
        xyzData = []
        for line in outData:
            if line[0] == '|':
                xyzData.append(line)
        # Remove the header lines
        while '*' not in xyzData[0]:
            xyzData.pop(0)
        xyzData.pop(0)
        # Remove final lines and ORCA formatting
        popList = []
        for lineNo, line in enumerate(xyzData):
            if ('*' in line) or ('end' in line) or (len(line.split('> ')[-1][:-1]) < 3):
                popList.append(lineNo)
            else:
                xyzData[lineNo] = line.split('> ')[-1]
        for lineNo in sorted(popList, reverse=True):
            xyzData.pop(lineNo)
        xyzData.insert(0, str(len(xyzData))+'\n')
        xyzData.insert(1, 'CommentLine\n')
        with open('./'+fileName.replace('.out', '.xyz'), 'w+') as xyzFile:
            xyzFile.writelines(xyzData)
