from hoomd_script import *

import numpy as np
import os
import sys
sys.path.append('../../code')
import helperFunctions


def multiHarmonicTorsion(theta, V0, V1, V2, V3, V4):
    # Definition of multiharmonic dihedral equation based on 5 input parameters to be used by HOOMD
    V = V0 + V1 * np.cos(theta) + V2 * ((np.cos(theta))**2) + V3 * ((np.cos(theta))**3) + V4 * ((np.cos(theta))**4)
    F = V1 * np.sin(theta) + 2 * V2 * np.cos(theta) * np.sin(theta) + 3 * V3 * ((np.cos(theta))**2) * np.sin(theta) + 4 * V4 * ((np.cos(theta))**3) * np.sin(theta)
    return (V, F)


def getFFCoeffs(FFName):
    masterFF = {}
    FFList = helperFunctions.loadFFXML(FFName)
    for FFType in list(FFList.keys()):
        if FFType not in list(masterFF.keys()):
            masterFF[FFType] = FFList[FFType]
        else:
            masterFF[FFType] += FFList[FFType]
    # finally, assign the expected variables to each value in the masterFF
    ljCoeffs = masterFF['lj']
    dpdCoeffs = masterFF['dpd']
    bondCoeffs = masterFF['bond']
    angleCoeffs = masterFF['angle']
    dihedralCoeffs = masterFF['dihedral']
    improperCoeffs = masterFF['improper']

    # set pair Coeffs
    lj = pair.lj(r_cut=2.5)
    #lj = hoomd.md.pair.lj(r_cut=2.5, nlist=nl)
    for index1, ljFF1 in enumerate(ljCoeffs):
        for index2, ljFF2 in enumerate(ljCoeffs):
            if index1 < index2:
                continue
            lj.pair_coeff.set(ljFF1[0], ljFF2[0], epsilon=np.sqrt(ljFF1[1] * ljFF2[1]), sigma=np.sqrt(ljFF1[2] * ljFF2[2]))

    # set bond coeffs
    harmonic_bond = bond.harmonic()
    for bondCoeff in bondCoeffs:
        # [k] = kcal mol^{-1} \aa^{-2} * episilon/sigma^{2}, [r0] = \aa * sigma^{2}
        harmonic_bond.bond_coeff.set(bondCoeff[0], k=bondCoeff[1], r0=bondCoeff[2])

    # set angle coeffs
    harmonic_angle = angle.harmonic()
    for angleCoeff in angleCoeffs:
    # [k] = kcal mol^{-1} rad^{-2} * epsilon, [t] = rad
        harmonic_angle.set_coeff(angleCoeff[0], k=angleCoeff[1], t0=angleCoeff[2])

    # set dihedral coeffs
    if len(dihedralCoeffs[0]) == 6:
        # Multiharmonic
        table_dihedral = dihedral.table(width=1000)
        for dihedralCoeff in dihedralCoeffs:
            table_dihedral.dihedral_coeff.set(dihedralCoeff[0], func=multiHarmonicTorsion, coeff=dict(V0=dihedralCoeff[1], V1=dihedralCoeff[2], V2=dihedralCoeff[3], V3=dihedralCoeff[4], V4=dihedralCoeff[5]))
    elif len(dihedralCoeffs[0]) == 5:
        # OPLS
        opls_dihedral = dihedral.opls()
        for dihedralCoeff in dihedralCoeffs:
            opls_dihedral.set_coeff(dihedralCoeff[0], k1 = dihedralCoeff[1], k2 = dihedralCoeff[2], k3 = dihedralCoeff[3], k4 = dihedralCoeff[4])
    else:
        raise SystemError("UNKNOWN DIHEDRAL TYPE")

    # set Improper Coeffs
    if len(improperCoeffs) > 0:
        harmonic_improper = improper.harmonic()
        for improperCoeff in improperCoeffs:
            harmonic_improper.improper_coeff.set(improperCoeff[0], k=improperCoeff[1], chi=improperCoeff[2])


if __name__ == "__main__":
    try:
        fileName = sys.argv[1]
        FFFileName = sys.argv[2]
        temperature = float(sys.argv[3])
        targetBoxSize = float(sys.argv[4])
    except:
        print("Please run this script using the following syntax:")
        print("python shrinkHOOMD1AA.py <MORPHOLOGY_FILE_LOC> <FORCEFIELD_FILE_LOC> <TEMPERATURE> <TARGET_BOX_SIZE>")
        exit()

    print("RUNNING", fileName)

    system = init.read_xml(filename=fileName)
    getFFCoeffs(FFFileName)

    allGroup = group.all()
    rigidGroup = group.rigid()
    nonRigidGroup = group.nonrigid()

    # Get the initial temperature of the simulation
    #hyphenLocs = helperFunctions.findIndex(fileName, '-')
    #temperature = fileName[hyphenLocs[3]+1:hyphenLocs[4]][1:]  # HARD CODED for the standard Jankowski naming nomenclature

    integrate.mode_standard(dt=0.001);
    rigidIntegrator = integrate.nvt_rigid(group=rigidGroup, T=temperature, tau=1.0)
    nonRigidIntegrator = integrate.nvt(group=nonRigidGroup, T=temperature, tau=1.0)

    run_time = 1e7

    dump.dcd(filename=fileName.replace(".xml", ".dcd"), overwrite=True, period=int(run_time/500))
    analyze.log(filename=fileName.replace(".xml", ".log"), quantities=['potential_energy'],
                period=int(run_time/1000), header_prefix='#', overwrite=True)

    # Get the initial box size dynamically
    initialMorphology = helperFunctions.loadMorphologyXML(fileName)
    update.box_resize(L = variant.linear_interp([(0, initialMorphology['lx']), (run_time, targetBoxSize)]))
    #update.box_resize(L = variant.linear_interp([(0, initialMorphology['lx']), (run_time, 85.18963051)]))  # HARD CODED for the ordered P3HT morphology volume size (atomistic, no scaling)
    run(run_time)
    dump.xml(filename="postshrink_" + fileName, all=True)
    run(run_time)
    dump.xml(filename="posteql_" + fileName, all=True)
