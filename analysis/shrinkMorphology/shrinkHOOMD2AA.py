import hoomd
import hoomd.md
import hoomd.deprecated

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
    nl = hoomd.md.nlist.cell()

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
    lj = hoomd.md.pair.lj(r_cut=2.5, nlist=nl)
    for index1, ljFF1 in enumerate(ljCoeffs):
        for index2, ljFF2 in enumerate(ljCoeffs):
            if index1 < index2:
                continue
            lj.pair_coeff.set(ljFF1[0], ljFF2[0], epsilon=np.sqrt(ljFF1[1] * ljFF2[1]), sigma=np.sqrt(ljFF1[2] * ljFF2[2]))

    # set bond coeffs
    harmonic_bond = hoomd.md.bond.harmonic()
    for bondCoeff in bondCoeffs:
        # [k] = kcal mol^{-1} \aa^{-2} * episilon/sigma^{2}, [r0] = \aa * sigma^{2}
        harmonic_bond.bond_coeff.set(bondCoeff[0], k=bondCoeff[1], r0=bondCoeff[2])

    # set angle coeffs
    harmonic_angle = hoomd.md.angle.harmonic()
    for angleCoeff in angleCoeffs:
    # [k] = kcal mol^{-1} rad^{-2} * epsilon, [t] = rad
        harmonic_angle.angle_coeff.set(angleCoeff[0], k=angleCoeff[1], t0=angleCoeff[2])

    # set dihedral coeffs
    table_dihedral = hoomd.md.dihedral.table(width=1000)
    for dihedralCoeff in dihedralCoeffs:
        table_dihedral.dihedral_coeff.set(dihedralCoeff[0], func=multiHarmonicTorsion, coeff=dict(V0=dihedralCoeff[1], V1=dihedralCoeff[2], V2=dihedralCoeff[3], V3=dihedralCoeff[4], V4=dihedralCoeff[5]))

    # set Improper Coeffs
    if len(improperCoeffs) > 0:
        harmonic_improper = hoomd.md.improper.harmonic()
        for improperCoeff in improperCoeffs:
            harmonic_improper.improper_coeff.set(improperCoeff[0], k=improperCoeff[1], chi=improperCoeff[2])


if __name__ == "__main__":
    FFFileName = 'FFP3HT.xml'
    filesToRun = []
    for fileName in os.listdir():
        if ('.xml' in fileName) and ('post' not in fileName) and ('FF' not in fileName):
            filesToRun.append(fileName)

    for fileName in filesToRun:
        hoomd.context.initialize("")
        system = hoomd.deprecated.init.read_xml(filename=fileName)
        getFFCoeffs(FFFileName)

        groupAll = hoomd.group.all()

        # Get the initial temperature of the simulation
        hyphenLocs = helperFunctions.findIndex(fileName, '-')
        temperature = fileName[hyphenLocs[3]+1:hyphenLocs[4]][1:]  # HARD CODED for the standard Jankowski naming nomenclature

        hoomd.md.integrate.mode_standard(dt=0.0001);
        integrator = hoomd.md.integrate.nvt(group=all, tau=1.0, kT=temperature)

        run_time = 1e7

        hoomd.dump.dcd(filename=fileName + ".dcd", period=int(run_time/500), overwrite=True)
        hoomd.analyze.log(filename=fileName + '.log', quantities=['potential_energy'],
                            period=int(run_time/1000), header_prefix='#', overwrite=True)

        # Get the initial box size dynamically
        initialMorphology = helperFunctions.loadMorphologyXML(fileName)
        hoomd.update.box_resize(L = hoomd.variant.linear_interp([(0, initialMorphology['lx']), (run_time, 28.39654350281)]))  # HARD CODED for the ordered P3HT morphology volume size
        hoomd.run(run_time)
        hoomd.deprecated.dump.xml(group=all, filename="postshrink_" + fileName, all=True)
        hoomd.run(run_time)
        hoomd.deprecated.dump.xml(group=all, filename="posteql_" + fileName, all=True)
