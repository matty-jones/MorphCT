import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import csv

from scipy.optimize import curve_fit


def loadData(CSVDir):
    hopTimeData = {}
    for fileName in os.listdir(CSVDir):
        if '.csv' in fileName:
            with open(CSVDir+'/'+fileName, 'r') as csvFile:
                csvReader = csv.reader(csvFile, delimiter=',')
                temp = fileName[:5]
                hopTimeData[temp] = []
                for row in csvReader:
                    hopTimeData[temp].append(float(row[0]))
    return hopTimeData




def plotHist(hopTimeData):
    for key in hopTimeData.keys():
        plt.hist(hopTimeData[key], bins=np.logspace(1E-21, 1E-10, 12))
        plt.gca().set_xscale('log')
        plt.xlabel('Hop Times (s)')
        plt.ylabel('Frequency')
        plt.savefig('./hist_'+key+'.pdf')
        plt.clf()
        print "Figure saved as", './hist_'+key+'.pdf'


if __name__ == "__main__":
    hopTimeData = loadData('./')
    plotHist(hopTimeData)
