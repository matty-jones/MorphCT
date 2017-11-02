import sys
import os

if __name__ == "__main__":
    files = os.listdir('./')
    for fileName in files:
        if '.inp' not in fileName:
            continue
        print(fileName)
        with open('./'+fileName, 'r') as inpFile:
            inData = inpFile.readlines()
        xyzData = []
        for line in inData:
            if line[0] == ' ':
                xyzData.append(line)
        xyzData.insert(0, str(len(xyzData))+'\n')
        xyzData.insert(1, 'CommentLine\n')
        with open('./'+fileName.replace('.inp', '.xyz'), 'w+') as xyzFile:
            xyzFile.writelines(xyzData)
