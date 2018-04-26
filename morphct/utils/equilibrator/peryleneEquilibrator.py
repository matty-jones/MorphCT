import sys
sys.path.append('../../code')
import helperFunctions
import copy
import argparse
import numpy as np
import hoomd
import hoomd.md
import hoomd.deprecated


def setCoeffs(morphologyDict):
    # Perylene: Sigma = 3.8 Ang, Epsilon = 0.1217 kCal/mol
    morphologyBonds = sorted(list(set([x[0] for x in morphologyDict['bond']])))
    morphologyAngles = sorted(list(set([x[0] for x in morphologyDict['angle']])))
    morphologyDihedrals = sorted(list(set([x[0] for x in morphologyDict['dihedral']])))
    morphologyImpropers = sorted(list(set([x[0] for x in morphologyDict['improper']])))


    ljnl = hoomd.md.nlist.cell()
    lj = hoomd.md.pair.lj(r_cut=2.5, nlist=ljnl)
    carbonTypes = ['C' + str(x + 1) for x in range(20)]
    for carbon1 in carbonTypes:
        for carbon2 in carbonTypes:
            lj.pair_coeff.set(carbon1, carbon2, epsilon=0.32, sigma=1.0)

    pppmnl = hoomd.md.nlist.cell()
    pppm = hoomd.md.charge.pppm(group=hoomd.group.charged(), nlist = pppmnl)
    pppm.set_params(Nx=64,Ny=64,Nz=64,order=6,rcut=2.70)

    harmonic_bond = hoomd.md.bond.harmonic()
    for bond in morphologyBonds:
        harmonic_bond.bond_coeff.set(bond, k=30000.0, r0=0.4)

    harmonic_angle = hoomd.md.angle.harmonic()
    for angle in morphologyAngles:
        harmonic_angle.angle_coeff.set(angle, k=380.0, t0=2.09)

    harmonic_dihedral = hoomd.md.dihedral.opls()
    for dihedral in morphologyDihedrals:
        harmonic_dihedral.dihedral_coeff.set(dihedral, k1=0.0, k2=50.0, k3=0.0, k4=0.0)


def modifyMorphology(morphologyDict, charges):
    print("Resetting Rigid Bodies...")
    morphologyDict['body'] = [-1] * len(morphologyDict['type'])
    if len(morphologyDict['bond']) == 0:
        print("Adding Constraints to System...")
        flexibleMorphologyDict = addConstraints(morphologyDict)
    else:
        flexibleMorphologyDict = copy.deepcopy(morphologyDict)
    print("Adding Charges to System...")
    chargedMorphologyDict = addCharges(flexibleMorphologyDict, charges)
    return chargedMorphologyDict


def addConstraints(morphologyDict):
    templateMorph = helperFunctions.loadMorphologyXML('./PEtemplate.xml')
    # One molecule of Perylene contains 20 atoms.
    # If we assume the molecules are coherent, we can add 20 to the constraint atomIDs for each molecule
    constraintTypes = ['bond', 'angle', 'dihedral', 'improper']
    for IDModifier in range(0, len(morphologyDict['type']), 20):
        # Iterate over constraint types
        for constraintType in constraintTypes:
            # Iterate over constraints
            for constraint in templateMorph[constraintType]:
                # Append the morphology with the new constraint, where the atom IDs have been modified
                morphologyDict[constraintType].append([constraint[0]] + [x + IDModifier for x in constraint[1:]])
    return morphologyDict


def addCharges(morphologyDict, charges):
    # Add the charges (zeroed out so that we can incrementally update them)
    morphologyDict['charge'] = [0 for _ in range(len(morphologyDict['type']))]
    for atomIndex, atomType in enumerate(morphologyDict['type']):
        morphologyDict['charge'][atomIndex] = charges[atomType]
    return morphologyDict


def initializeVelocities(snapshot, temperature):
    v = np.random.random((len(snapshot.particles.velocity), 3))
    v -= 0.5
    meanv = np.mean(v, 0)
    meanv2 = np.mean(v ** 2, 0)
    fs = np.sqrt(temperature / meanv2)
    # Shift the velocities such that the average is zero
    v = (v - meanv)
    # Scale the velocities to match the required temperature
    v *= fs
    # Assign the velocities for this MD phase
    snapshot.particles.velocity[:] = v[:]
    return snapshot


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--temperature", type=float, required=True, help="The temperature of the system")
    #parser.add_argument("-p", "--perylothiophene", action="store_true", required=False, help="If present, allow periodic connections to add chromophores to stacks, as well as non-periodic connections (this usually reduces the number of stacks in the system). Defaults to False.")
    args, fileList = parser.parse_known_args()

    # OLD CHARGES
    #charges = {'CP': 3.880, 'CN': -5.820}
    # PE Charges given by PECharges_NWChem.out in this directory
    # C Atoms (+ bonded Hydrogens):
    # 1 + 21 = -0.190034 + 0.131579 = -0.058455
    # 2 = 0.151918
    # 3 + 22 = -0.191069 + 0.150689 = -0.040380
    # 4 + 23 = -0.124397 + 0.124275 = -0.000122
    # 5 = 0.057477
    # 6 = -0.043335
    # 7 + 24 = -0.198055 + 0.146827 = -0.051228
    # 8 + 25 = -0.159231 + 0.138569 = -0.020662
    # 9 = 0.052654
    # 10 + 26 = -0.170768 + 0.143989 = -0.026779
    # 11 = -0.010152
    # 12 = 0.038679
    # 13 = -0.003800
    # 14 + 27 = -0.149203 + 0.133416 = -0.015787
    # 15 + 28 = -0.120864 + 0.121221 = 0.000357
    # 16 + 29 = -0.244458 + 0.164674 = -0.079784
    # 17 = 0.145754
    # 18 + 30 = -0.232718 + 0.151866 = -0.080852
    # 19 + 31 = -0.085985 + 0.106369 = 0.020384
    # 20 + 32 = -0.180134 + 0.144250 = -0.035884
    chargesE = {'C1': -0.058455, 'C2': 0.151915, 'C3': -0.040380,
                'C4': -0.000122, 'C5': 0.057477, 'C6': -0.043335,
                'C7': -0.051228, 'C8': -0.020662, 'C9': 0.052654,
                'C10': -0.026779, 'C11': -0.010152, 'C12': 0.038679,
                'C13': -0.003800, 'C14': -0.015787, 'C15': 0.000357,
                'C16': -0.079784, 'C17': 0.145754, 'C18': -0.080852,
                'C19': 0.020384, 'C20': -0.035884}
    # Note: I subtracted 3E-6 to the 'C2' charge to make it sum to zero (rounding error in NWChem output)
    charges = {}
    # Convert DFT outputs to HOOMD dimensionless units. Epsilon = 0.122 kCal/mol =
    # 0.510 kJ/mol = 8.47920266E-22 J
    # Sigma = 3.8E-10 m
    # E0 = 8.85E-12 [SI]
    e0 = 8.85418782E-12
    sigma = 3.8E-10
    epsilon = 8.47920266E-22
    chargeFactor = 1.60217662E-19 / (4 * np.pi * e0 * sigma * epsilon)**0.5
    for atomType, charge in chargesE.items():
        charges[atomType] = charge * chargeFactor
    for atomType, charge in charges.items():
        print(atomType, charge)
    #print("Average Positive charge =", np.average(np.array([val for val in charges.values() if val > 0.0])))
    #print("Average Negative charge =", np.average(np.array([val for val in charges.values() if val < 0.0])))
    #exit()
    # Masses normalized for M_{sulfur}, as in PT.
    masses = {}
    for DFTCHIndex in [1, 3, 4, 7, 8, 10, 14, 15, 16, 18, 19, 20]:
        masses['C' + str(DFTCHIndex)] = 0.406
    for DFTCIndex in [2, 5, 6, 9, 11, 12, 13, 17]:
        masses['C' + str(DFTCIndex)] = 0.375
    chargeIncrements = 20
    chargeTimesteps = 10000
    run_time = 3e8

    for fileName in fileList:
        print("Running", fileName)
        morphologyDict = helperFunctions.loadMorphologyXML(fileName)
        morphologyDict = modifyMorphology(morphologyDict, charges)
        newFileName = fileName.split('.')[0] + '_wConsCharge.xml'
        helperFunctions.writeMorphologyXML(morphologyDict, newFileName)

        hoomd.context.initialize("")
        system = hoomd.deprecated.init.read_xml(filename=newFileName)
        # Set the correct atom Masses
        for atom in system.particles:
            atom.mass = masses[atom.type]

        snapshot = system.take_snapshot()
        # Assign the required velocities based on the requested temperature
        updatedSnapshot = initializeVelocities(snapshot, args.temperature)
        # Finally, restore the snapshot
        system.restore_snapshot(updatedSnapshot)

        setCoeffs(morphologyDict)

        hoomd.md.integrate.mode_standard(dt=0.001);
        integrator = hoomd.md.integrate.nvt(group=hoomd.group.all(), tau=1.0, kT=args.temperature)
        #hoomd.md.nlist.reset_exclusions(exclusions = ['bond', 'angle', 'dihedral', 'body'])

        hoomd.dump.dcd(filename=fileName.split('.')[0] + ".dcd", period=int(run_time/500), overwrite=True)
        logQuantities = ['temperature', 'pressure', 'volume', 'potential_energy', 'kinetic_energy', 'pair_lj_energy', 'pair_ewald_energy', 'pppm_energy', 'bond_harmonic_energy', 'angle_harmonic_energy', 'dihedral_opls_energy']
        hoomd.analyze.log(filename=fileName.split('.')[0] + ".log", quantities=logQuantities,
                            period=int(run_time/10000), header_prefix='#', overwrite=True)
        # Now incrementally ramp the charges
        for chargePhase in range(chargeIncrements + 1):
            print("Incrementing charge phase", chargePhase, "of", chargeIncrements + 1)
            for atom in system.particles:
                oldCharge = copy.deepcopy(atom.charge)
                atom.charge = charges[atom.type] * (chargePhase / float(chargeIncrements))
            hoomd.run(chargeTimesteps)

        # Get the initial box size dynamically
        hoomd.run_upto(run_time)
        hoomd.deprecated.dump.xml(group=hoomd.group.all(), filename="relaxed_" + fileName, all=True)
