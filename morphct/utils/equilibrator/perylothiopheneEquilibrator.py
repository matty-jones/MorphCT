import sys
sys.path.append('../../code')
import helperFunctions
import copy
import argparse
import hoomd
import hoomd.md
import hoomd.deprecated
import numpy as np


def setCoeffs(morphologyDict):
    # Perylothiophene: Sigma = 3.8 Ang, Epsilon = 0.358 kCal/mol
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
        lj.pair_coeff.set(carbon1, 'S', epsilon=0.57, sigma=0.96)
    lj.pair_coeff.set('S', 'S', epsilon=1.0, sigma=0.92)

    pppmnl = hoomd.md.nlist.cell()
    pppm = hoomd.md.charge.pppm(group=hoomd.group.charged(), nlist = pppmnl)
    pppm.set_params(Nx=64,Ny=64,Nz=64,order=6,rcut=2.70)

    harmonic_bond = hoomd.md.bond.harmonic()
    for bond in morphologyBonds:
        if bond.count('C') == 1:
            harmonic_bond.bond_coeff.set(bond, k=30000.0, r0=0.45)
        elif bond.count('C') == 2:
            harmonic_bond.bond_coeff.set(bond, k=30000.0, r0=0.4)

    harmonic_angle = hoomd.md.angle.harmonic()
    for angle in morphologyAngles:
        if angle.count('C') == 2:
            harmonic_angle.angle_coeff.set(angle, k=380.0, t0=1.60)
        elif angle.count('C') == 3:
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
    templateMorph = helperFunctions.loadMorphologyXML('./relabeled_PTTemplate.xml')
    # One molecule of Perylothiophene contains 21 atoms.
    # If we assume the molecules are coherent, we can add 21 to the constraint atomIDs for each molecule
    constraintTypes = ['bond', 'angle', 'dihedral', 'improper']
    for IDModifier in range(0, len(morphologyDict['type']), 21):
        # Iterate over constraint types
        for constraintType in constraintTypes:
            # Iterate over constraints
            for constraint in templateMorph[constraintType]:
                # Append the morphology with the new constraint, where the atom IDs have been modified
                morphologyDict[constraintType].append([constraint[0]] + [x + IDModifier for x in constraint[1:]])
    return morphologyDict


def addCharges(morphologyDict, charges):
    # Firstly, need to update the atom types based on the number of bonds they have.
    fixedAtomTypesDict = fixAtomTypes(morphologyDict)
    # Then add the charges (zeroed out so that we can incrementally update them)
    fixedAtomTypesDict['charge'] = [0 for _ in range(len(fixedAtomTypesDict['type']))]
    for atomIndex, atomType in enumerate(fixedAtomTypesDict['type']):
        fixedAtomTypesDict['charge'][atomIndex] = charges[atomType]
    #print(sorted(list(set([x[0] for x in fixedAtomTypesDict['bond']]))))
    #print(sorted(list(set([x[0] for x in fixedAtomTypesDict['angle']]))))
    #print(sorted(list(set([x[0] for x in fixedAtomTypesDict['dihedral']]))))
    #print(sorted(list(set([x[0] for x in fixedAtomTypesDict['improper']]))))
    #exit()
    return fixedAtomTypesDict


def fixAtomTypes(morphologyDict):
    # Iterate through the atom typelist to give the correct atom types
    print(len(morphologyDict['type']))
    for AAID, typeName in enumerate(morphologyDict['type']):
        if typeName == 'C':
            morphologyDict['type'][AAID] = 'C' + str(AAID%21 + 1)
    # Finally, fix all the constraint names
    constraintTypes = ['bond', 'angle', 'dihedral', 'improper']
    for constraintType in constraintTypes:
        for constraintIndex, constraint in enumerate(morphologyDict[constraintType]):
            newConstraint = copy.deepcopy(constraint)
            newConstraint[0] = '-'.join([morphologyDict['type'][atomID] for atomID in newConstraint[1:]])
            morphologyDict[constraintType][constraintIndex] = newConstraint
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
    # PT Charges given by PTCharges_NWChem.out in this directory
    # Calculated by: CH = <CH>, CC = <C>, S = S, where <> denotes average and CH is the carbon + bonded hydrogen for relevant atoms
    # C Atoms (+ bonded Hydrogens):
    # 1 + 22 = -0.232266 + 0.149556 = -0.082710
    # 2 = 0.157736
    # 3 + 23 = -0.200920 + 0.160898 = -0.040022
    # 4 + 24 = -0.142969 + 0.132762 = -0.010207
    # 5 = 0.016040
    # 6 = 0.069730
    # 7 + 25 = -0.238319 + 0.151989 = -0.086330
    # 8 + 26 = -0.234022 + 0.198674 = -0.035348
    # 9 = -0.111277
    # 10 = 0.152129
    # 11 = 0.055825
    # 12 = 0.086899
    # 13 = -0.126481
    # 14 + 27 = -0.240882 + 0.181823 = -0.059059
    # 15 + 28 = -0.155182 + 0.140634 = -0.014548
    # 16 + 29 = -0.177484 + 0.142082 = -0.035402
    # 17 = 0.112991
    # 18 + 30 = -0.255068 + 0.152194 = -0.102874
    # 19 + 31 = -0.228936 + 0.184658 = -0.044278
    # 20 = 0.182532
    # S Atoms:
    # 21 = -0.085347
    chargesE = {'C1': -0.082710, 'C2': 0.157736, 'C3': -0.040022,
                'C4': -0.010207, 'C5': 0.016040, 'C6': 0.069730,
                'C7': -0.086330, 'C8': -0.035348, 'C9': -0.111277,
                'C10': 0.152129, 'C11': 0.055825, 'C12': 0.086899,
                'C13': -0.126481, 'C14': -0.059059, 'C15': -0.014548,
                'C16': -0.035402, 'C17': 0.112991, 'C18': -0.102874,
                'C19': -0.044278, 'C20': 0.182532, 'S': -0.085346}
    # Note: I added 1E-6 to the 'S' charge to make it sum to zero (rounding error in NWChem output)
    charges = {}
    # Convert DFT outputs to HOOMD dimensionless units. Epsilon = 0.358 kCal/mol = 
    # 1.498 kJ/mol = 2.48727492E-21 J
    # Sigma = 3.8E-10 m
    # E0 = 8.85E-12 [SI]
    e0 = 8.85418782E-12
    sigma = 3.8E-10
    epsilon = 2.48727492E-21
    chargeFactor = 1.60217662E-19 / (4 * np.pi * e0 * sigma * epsilon)**0.5
    for atomType, charge in chargesE.items():
        charges[atomType] = charge * chargeFactor
    masses = {'S': 1.000}
    for DFTCHIndex in [1, 3, 4, 7, 8, 14, 16, 18, 19]:
        masses['C' + str(DFTCHIndex)] = 0.406
    for DFTCIndex in [2, 5, 6, 9, 10, 11, 12, 13, 15, 17, 20]:
        masses['C' + str(DFTCIndex)] = 0.375
    chargeIncrements = 20
    chargeTimesteps = 10000
    run_time = 3e8

    for fileName in fileList:
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
