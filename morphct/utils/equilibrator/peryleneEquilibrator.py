import sys
sys.path.append('../../code')
import helperFunctions
import copy
import argparse
import numpy as np
import hoomd
import hoomd.md
import hoomd.deprecated


def setCoeffs():
    # Perylene: Sigma = 3.8 Ang, Epsilon = 0.1217 kCal/mol
    ljnl = hoomd.md.nlist.cell()
    lj = hoomd.md.pair.lj(r_cut=2.5, nlist=ljnl)
    lj.pair_coeff.set('CP', 'CP', epsilon=1.0, sigma=1.0)
    lj.pair_coeff.set('CN', 'CN', epsilon=1.0, sigma=1.0)
    lj.pair_coeff.set('CP', 'CN', epsilon=1.0, sigma=1.0)

    pppmnl = hoomd.md.nlist.cell()
    pppm = hoomd.md.charge.pppm(group=hoomd.group.charged(), nlist = pppmnl)
    pppm.set_params(Nx=64,Ny=64,Nz=64,order=6,rcut=2.70)

    harmonic_bond = hoomd.md.bond.harmonic()
    harmonic_bond.bond_coeff.set('CP-CP', k=30000.0, r0=0.4)
    harmonic_bond.bond_coeff.set('CP-CN', k=30000.0, r0=0.4)
    harmonic_bond.bond_coeff.set('CN-CP', k=30000.0, r0=0.4)
    harmonic_bond.bond_coeff.set('CN-CN', k=30000.0, r0=0.4)

    harmonic_angle = hoomd.md.angle.harmonic()
    harmonic_angle.angle_coeff.set('CN-CN-CN', k=380.0, t0=2.09)
    harmonic_angle.angle_coeff.set('CN-CN-CP', k=380.0, t0=2.09)
    harmonic_angle.angle_coeff.set('CN-CP-CP', k=380.0, t0=2.09)
    harmonic_angle.angle_coeff.set('CP-CN-CN', k=380.0, t0=2.09)
    harmonic_angle.angle_coeff.set('CP-CN-CP', k=380.0, t0=2.09)
    harmonic_angle.angle_coeff.set('CP-CP-CN', k=380.0, t0=2.09)
    harmonic_angle.angle_coeff.set('CP-CP-CP', k=380.0, t0=2.09)

    harmonic_dihedral = hoomd.md.dihedral.opls()
    harmonic_dihedral.dihedral_coeff.set('CN-CN-CN-CN', k1=0.0, k2=50.0, k3=0.0, k4=0.0)
    harmonic_dihedral.dihedral_coeff.set('CN-CN-CN-CP', k1=0.0, k2=50.0, k3=0.0, k4=0.0)
    harmonic_dihedral.dihedral_coeff.set('CN-CN-CP-CP', k1=0.0, k2=50.0, k3=0.0, k4=0.0)
    harmonic_dihedral.dihedral_coeff.set('CN-CP-CP-CP', k1=0.0, k2=50.0, k3=0.0, k4=0.0)
    harmonic_dihedral.dihedral_coeff.set('CP-CN-CN-CN', k1=0.0, k2=50.0, k3=0.0, k4=0.0)
    harmonic_dihedral.dihedral_coeff.set('CP-CN-CN-CP', k1=0.0, k2=50.0, k3=0.0, k4=0.0)
    harmonic_dihedral.dihedral_coeff.set('CP-CN-CP-CP', k1=0.0, k2=50.0, k3=0.0, k4=0.0)
    harmonic_dihedral.dihedral_coeff.set('CP-CP-CN-CN', k1=0.0, k2=50.0, k3=0.0, k4=0.0)
    harmonic_dihedral.dihedral_coeff.set('CP-CP-CN-CP', k1=0.0, k2=50.0, k3=0.0, k4=0.0)
    harmonic_dihedral.dihedral_coeff.set('CP-CP-CP-CN', k1=0.0, k2=50.0, k3=0.0, k4=0.0)


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
    # Firstly, need to update the atom types based on the number of bonds they have.
    fixedAtomTypesDict = fixAtomTypes(morphologyDict)
    # Then add the charges (zeroed out so that we can incrementally update them)
    fixedAtomTypesDict['charge'] = [0 for _ in range(len(fixedAtomTypesDict['type']))]
    for atomIndex, atomType in enumerate(fixedAtomTypesDict['type']):
        if atomType == 'CP':
            fixedAtomTypesDict['charge'][atomIndex] = charges[atomType]
        elif atomType == 'CN':
            fixedAtomTypesDict['charge'][atomIndex] = charges[atomType]
    return fixedAtomTypesDict


def fixAtomTypes(morphologyDict):
    # If the carbon has 2 bonds, it is an `outer' atom that should have a positive charge ('CP')
    # If the carbon has 3 bonds, it is an `inner' atom that should have a negative charge ('CN')
    bondNumbers = {}
    for bond in morphologyDict['bond']:
        for atomID in bond[1:]:
            if atomID not in bondNumbers:
                bondNumbers[atomID] = 1
            else:
                bondNumbers[atomID] += 1
    positiveIDs = [atomID for atomID, bondNumber in bondNumbers.items() if bondNumber == 2]
    negativeIDs = [atomID for atomID, bondNumber in bondNumbers.items() if bondNumber == 3]
    for index in positiveIDs:
        morphologyDict['type'][index] = 'CP'
    for index in negativeIDs:
        morphologyDict['type'][index] = 'CN'
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

    charges = {'CP': 3.880, 'CN': -5.820}
    masses = {'CP': 1.000, 'CN': 0.923}
    chargeIncrements = 20
    chargeTimesteps = 10000
    run_time = 1e7

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

        setCoeffs()

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
