import pickle
import os
import subprocess as sp
import sys
from scipy.sparse import lil_matrix
import glob
import re
from fileinput import input
from shutil import copyfile

def create_results():
    sp.call("ls -ltr > tmp.txt", shell=True)
    cores_list = []
    for core in glob.glob("KMCData_*.pickle"):
        cores_list.append(re.findall("KMCData_(.*).pickle", core)[0])
    keep_list = []
    for core in cores_list:
        select_list = []
        with open("tmp.txt", 'r') as f:
            for line in f:
                if 'KMCslot' in line and 'Results_'+core in line:
                    select_list.append([line.split()[-1]])
        keep_list.append(select_list[-1])
    for keeper in zip(cores_list, keep_list):
        new_name = "KMCResults_"+str(keeper[0])+".pickle"
        copyfile(str(keeper[1][0]), new_name)
    os.remove('tmp.txt')

if __name__ == "__main__":
    create_results()
    outputDir = sys.argv[1]
    print("Combining outputs...")
    combinedData = {}
    pickleFiles = []
    for fileName in os.listdir(outputDir):
        if ("KMCResults_" in fileName):
            pickleFiles.append(outputDir + '/' + fileName)
    pickleFiles = sorted(pickleFiles)
    print("%d pickle files found to combine!" % (len(pickleFiles)))
    for fileName in pickleFiles:
        # The pickle was repeatedly dumped to, in order to save time.
        # Each dump stream is self-contained, so iteratively unpickle to add the new data.
        with open(fileName, 'rb') as pickleFile:
            pickledData = pickle.load(pickleFile)
            for key, val in pickledData.items():
                try:
                    if val is None:
                        continue
                    if key not in combinedData:
                        combinedData[key] = val
                    else:
                        combinedData[key] += val
                except AttributeError:
                    pass
    # Write out the combined data
    print("Writing out the combined pickle file...")
    with open(outputDir + '/KMCResults.pickle', 'wb+') as pickleFile:
        pickle.dump(combinedData, pickleFile)
    print("Complete data written to", outputDir + "/KMCResults.pickle.")
