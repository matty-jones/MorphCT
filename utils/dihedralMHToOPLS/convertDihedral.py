import numpy as np
import matplotlib.pyplot as plt


def MH(theta, V0, V1, V2, V3, V4):
    # Definition of multiharmonic dihedral equation based on 5 input parameters to be used by HOOMD
    V = V0 + (V1 * np.cos(theta)) + (V2 * ((np.cos(theta))**2)) + (V3 * ((np.cos(theta))**3)) + (V4 * ((np.cos(theta))**4))
    return V


def OPLS(theta, k1, k2, k3, k4):
    # Conventional OPLS fourier series
    V = 0.5 * ((k1 * (1 + np.cos(theta))) + (k2 * (1 - np.cos(2 * theta))) + (k3 * (1 + np.cos(3 * theta))) + (k4 * (1 - np.cos(4 * theta))))
    return V


def OPLStoMH(k1, k2, k3, k4):
    V0 = 0.5 * (k1 + (2 * k2) + k3)
    V1 = 0.5 * (k1 - (3 * k3))
    V2 = (4 * k4) - k2
    V3 = 2 * k3
    V4 = -(4 * k4)
    return V0, V1, V2, V3, V4


def MHtoOPLS(V0, V1, V2, V3, V4):
    k1 = (2 * V1) + (1.5 * V3)
    k2 = -V4 - V2
    k3 = 0.5 * V3
    k4 = -0.25 * V4
    return k1, k2, k3, k4


if __name__ == "__main__":
    plotting = True
    theta = np.arange(0, 2*np.pi, np.pi/32.)
    k1Vals = [0.0, 0.0, 0.0, 0.0, 0.0, 1.3, 0.0, 0.0, 4.669, 0.0, 0.0, 1.3, 0.0, 0.0, 0.0]
    k2Vals = [14.5, 14.5, 14.5, 14.5, 14.5, -0.05, 0.0, 0.0, 5.124, 5.124, 0.0, -0.05, 0.0, 14.5, 14.5]
    k3Vals = [0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.3, 0.3, 0.0, 0.0, 0.198, 0.2, 0.3, 0.0, 0.0]
    k4Vals = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    V0Vals = []
    V1Vals = []
    V2Vals = []
    V3Vals = []
    V4Vals = []

    for index, k1 in enumerate(k1Vals):
        k2 = k2Vals[index]
        k3 = k3Vals[index]
        k4 = k4Vals[index]
        actualPot = OPLS(theta, k1, k2, k3, k4)
        V0, V1, V2, V3, V4 = OPLStoMH(k1, k2, k3, k4)
        V0Vals.append(V0)
        V1Vals.append(V1)
        V2Vals.append(V2)
        V3Vals.append(V3)
        V4Vals.append(V4)

        if plotting is True:
            calcPot = MH(theta, V0, V1, V2, V3, V4)

            plt.figure()
            plt.plot(theta, actualPot, 'r', linewidth = 4.0)
            plt.plot(theta, calcPot, 'b', linewidth = 2.0)
            plt.savefig('./combined_forward_' + str(index) + '.pdf')

    k1Vals = []
    k2Vals = []
    k3Vals = []
    k4Vals = []

    for index, V0 in enumerate(V0Vals):
        V1 = V1Vals[index]
        V2 = V1Vals[index]
        V3 = V1Vals[index]
        V4 = V1Vals[index]
        actualPot = MH(theta, V0, V1, V2, V3, V4)
        k1, k2, k3, k4 = MHtoOPLS(V0, V1, V2, V3, V4)
        k1Vals.append(k1)
        k2Vals.append(k2)
        k3Vals.append(k3)
        k4Vals.append(k4)

        if plotting is True:
            calcPot = OPLS(theta, k1, k2, k3, k4)

            plt.figure()
            plt.plot(theta, actualPot, 'r', linewidth = 4.0)
            plt.plot(theta, calcPot, 'b', linewidth = 2.0)
            plt.savefig('./combined_backward_' + str(index) + '.pdf')



    #for index, k1 in enumerate(k1Vals):
    #    print "\n"
    #    print "The dihedral with OPLS Coeffs: %.4f, %.4f, %.4f, %.4f" % (k1, k2Vals[index], k3Vals[index], k4Vals[index])
    #    print "Should be replaced with these MultiHarmonic Coeffs: %.4f, %.4f, %.4f, %.4f, %.4f" % (V0Vals[index], V1Vals[index], V2Vals[index], V3Vals[index], V4Vals[index])
    #    raw_input("PAUSE...")

