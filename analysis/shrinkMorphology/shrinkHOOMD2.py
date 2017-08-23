import hoomd
import hoomd.md
import hoomd.deprecated

import os
import sys
sys.path.append('../../code')
import helperFunctions


def setCoeffs():
    nl = hoomd.md.nlist.cell()
    lj = hoomd.md.pair.lj(r_cut=2.5, nlist=nl)

    # Backbone-Backbone
    lj.pair_coeff.set('A', 'A', epsilon=2.0, sigma=1.0)

    # Backbone-Side Chain
    lj.pair_coeff.set('A', 'B', epsilon=0, sigma=1.0)
    lj.pair_coeff.set('A', 'C', epsilon=0, sigma=1.0)

    # Side Chain-Side Chain
    lj.pair_coeff.set('B', 'C', epsilon=0.1, sigma=1.0)
    lj.pair_coeff.set('C', 'C', epsilon=0.1, sigma=1.0)
    lj.pair_coeff.set('B', 'B', epsilon=0.1, sigma=1.0)

    harmonic_bond = hoomd.md.bond.harmonic()

    harmonic_bond.bond_coeff.set('bondA', k=100, r0=1.4)
    harmonic_bond.bond_coeff.set('bondB', k=100, r0=1.4)
    harmonic_bond.bond_coeff.set('bondC', k=100, r0=1.4)

    harmonic_angle = hoomd.md.angle.harmonic()

    harmonic_angle.angle_coeff.set('angleA', k=60.0, t0=3.14)
    harmonic_angle.angle_coeff.set('angleB', k=12.0, t0=2.13)
    harmonic_angle.angle_coeff.set('angleC', k=12.0, t0=2.13)
    harmonic_angle.angle_coeff.set('angleD', k=12.0, t0=3.14)

    harmonic_dihedral = hoomd.md.dihedral.harmonic()

    harmonic_dihedral.dihedral_coeff.set('dihedralA', k=45.0, d=1, n=1)


if __name__ == "__main__":
    filesToRun = []
    for fileName in os.listdir():
        if ('.xml' in fileName) and ('post' not in fileName):
            filesToRun.append(fileName)

    for fileName in filesToRun:
        hoomd.context.initialize("")
        system = hoomd.deprecated.init.read_xml(filename=fileName)
        all = hoomd.group.all()
        groupA = hoomd.group.type(name='a-particles', type='A')
        groupB = hoomd.group.type(name='b-particles', type='B')
        groupC = hoomd.group.type(name='c-particles', type='C')

        for p in groupB:
            p.mass = 0.518

        for p in groupC:
            p.mass = 0.530

        setCoeffs()  # HARD CODED for P3HT

        # Get the initial temperature of the simulation
        hyphenLocs = helperFunctions.findIndex(fileName, '-')
        temperature = fileName[hyphenLocs[3]+1:hyphenLocs[4]][1:]  # HARD CODED for the standard Jankowski naming nomenclature
        exit()

        hoomd.md.integrate.mode_standard(dt=0.001);
        integrator = hoomd.md.integrate.nvt(group=all, tau=1.0, kT=tmperature)

        run_time = 1e6

        hoomd.dump.dcd(filename="trajectory.dcd", period=int(run_time/100), overwrite=True)
        hoomd.analyze.log(filename='mylog.log', quantities=['potential_energy'],
                            period=int(run_time/1000), header_prefix='#', overwrite=True)

        # Get the initial box size dynamically
        initialMorphology = helperFunctions.loadMorphologyXML(fileName)
        hoomd.update.box_resize(L = hoomd.variant.linear_interp([(0, initialMorphology['lx']), (run_time, 28.39654350281)]))  # HARD CODED for the ordered P3HT morphology volume size
        hoomd.run(run_time)
        hoomd.deprecated.dump.xml(group=all, filename="postshrink_" + fileName, all=True)
        hoomd.run(run_time)
        hoomd.deprecated.dump.xml(group=all, filename="posteql_" + fileName, all=True)
