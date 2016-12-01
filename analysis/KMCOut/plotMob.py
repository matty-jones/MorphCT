import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import csv
import scipy.optimize
import scipy.stats


elementaryCharge = 1.60217657E-19  # C
kB = 1.3806488E-23  # m^{2} kg s^{-2} K^{-1}
temperature = 290  # K

if __name__ == "__main__":
    directory = sys.argv[1]

