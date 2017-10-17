import sys
sys.path.append('../../code')
import helperFunctions
import matplotlib.pyplot as plt
import numpy as np
try:
    import mpl_toolkits.mplot3d as p3
except ImportError:
    print("Could not import 3D plotting engine, calling the plotMolecule3D function will result in an error!")


def plotConnections(chromophoreList, simExtent, inputFile):
    # A complicated function that shows connections between carriers in 3D that carriers prefer to hop between.
    # Connections that are frequently used are highlighted in black, whereas rarely used connections are more white.
    # Find a good normalisation factor
    minimum = 9E99
    maximum = 0
    for chromophore in chromophoreList:
        for neighbourTI in chromophore.neighboursTI:
            if neighbourTI is None:
                neighbourTI = 0
            elif neighbourTI > maximum:
                maximum = neighbourTI
            elif neighbourTI < minimum:
                minimum = neighbourTI
    normalizeTo = float(maximum)
    # Try to get the colour map first
    plt.gcf()
    plt.clf()
    # Now for the actual plot
    fig = plt.gcf()
    ax = p3.Axes3D(fig)
    for chromo1 in chromophoreList:
        for index, chromo2Details in enumerate(chromo1.neighbours):
            coords1 = chromo1.posn
            coords2 = chromophoreList[chromo2Details[0]].posn
            value = chromo1.neighboursTI[index]
            if (not np.array_equal(chromo2Details[1], [0, 0, 0])) or (chromo1.ID > chromophoreList[chromo2Details[0]].ID):
                continue
            # Only plot connections between chromophores in the same image
            line = [coords2[0] - coords1[0], coords2[1] - coords1[1], coords2[2] - coords2[1]]
            if (np.abs(coords2[0] - coords1[0]) < simExtent[0] / 2.0) and (np.abs(coords2[1] - coords1[1]) < simExtent[1] / 2.0) and (np.abs(coords2[2] - coords1[2]) < simExtent[2] / 2.0):
                ax.plot([coords1[0], coords2[0]], [coords1[1], coords2[1]], [coords1[2], coords2[2]], c = 'k', linewidth = 0.5)
    fileName = '3d_' + inputFile[:-7] + '.pdf'
    plt.savefig('./' + fileName, bbox_inches='tight')
    print("Figure saved as", "./" + fileName)
    plt.clf()


if __name__ == "__main__":
    pickleFile = sys.argv[1]
    data = helperFunctions.loadPickle(pickleFile)
    simExtent = [data[0]['lx'], data[0]['ly'], data[0]['lz']]
    plotConnections(data[-1], simExtent, pickleFile)
