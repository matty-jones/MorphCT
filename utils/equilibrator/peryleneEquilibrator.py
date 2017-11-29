import sys
sys.path.append('../../code')
import helperFunctions
import copy
import argparse
import hoomd
import hoomd.md
import hoomd.deprecated


def setCoeffs():
    nl = hoomd.md.nlist.cell()
    lj = hoomd.md.pair.lj(r_cut=2.5, nlist=nl)
    lj.pair_coeff.set('CP', 'CP', epsilon=1.0, sigma=1.0)
    lj.pair_coeff.set('CN', 'CN', epsilon=1.0, sigma=1.0)
    lj.pair_coeff.set('CP', 'CN', epsilon=1.0, sigma=1.0)

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


def modifyMorphology(morphologyDict):
    if len(morphologyDict['bond']) == 0:
        print("Adding Constraints to System...")
        flexibleMorphologyDict = addConstraints(morphologyDict)
    else:
        flexibleMorphologyDict = copy.deepcopy(morphologyDict)
    print("Adding Charges to System...")
    chargedMorphologyDict = addCharges(flexibleMorphologyDict)
    return chargedMorphologyDict


def addConstraints(morphologyDict):
    templateMorph = helperFunctions.loadMorphologyXML('./template.xml')
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


def addCharges(morphologyDict):
    # Firstly, need to update the atom types based on the number of bonds they have.
    fixedAtomTypesDict = fixAtomTypes(morphologyDict)
    # Then add the charges
    # Zero out the old charges first
    fixedAtomTypesDict['charge'] = [0 for _ in range(len(fixedAtomTypesDict['type']))]
    for atomIndex, atomType in enumerate(fixedAtomTypesDict['type']):
        if atomType == 'CP':
            fixedAtomTypesDict['charge'][atomIndex] = 3.880
        elif atomType == 'CN':
            fixedAtomTypesDict['charge'][atomIndex] = -5.820
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--temperature", type=float, required=True, help="The temperature of the system")
    #parser.add_argument("-p", "--perylothiophene", action="store_true", required=False, help="If present, allow periodic connections to add chromophores to stacks, as well as non-periodic connections (this usually reduces the number of stacks in the system). Defaults to False.")
    args, fileList = parser.parse_known_args()

    for fileName in fileList:
        morphologyDict = helperFunctions.loadMorphologyXML(fileName)
        morphologyDict = modifyMorphology(morphologyDict)
        newFileName = fileName.split('.')[0] + '_wConsCharge.xml'
        helperFunctions.writeMorphologyXML(morphologyDict, newFileName)

        hoomd.context.initialize("")
        system = hoomd.deprecated.init.read_xml(filename=newFileName)
        all = hoomd.group.all()
        CPAtoms = hoomd.group.type(name='positive-carbons', type='CP')
        CNAtoms = hoomd.group.type(name='negative-carbons', type='CN')
        pppm = hoomd.md.charge.pppm(group=all)

        for p in CPAtoms:
            p.mass = 1.000

        for p in CNAtoms:
            p.mass = 0.923

        setCoeffs()

        hoomd.md.integrate.mode_standard(dt=0.001);
        integrator = hoomd.md.integrate.nvt(group=all, tau=1.0, kT=args.temperature)

        run_time = 1e7

        hoomd.dump.dcd(filename=fileName.split('.')[0] + ".dcd", period=int(run_time/500), overwrite=True)
        logQuantities = ['temperature', 'pressure', 'volume', 'potential_energy', 'kinetic_energy', 'pair_lj_energy', 'pair_ewald_energy', 'pppm_energy', 'bond_harmonic_energy', 'angle_harmonic_energy', 'dihedral_opls_energy']
        hoomd.analyze.log(filename=fileName.split('.')[0] + ".log", quantities=logQuantities,
                            period=int(run_time/1000), header_prefix='#', overwrite=True)

        # Get the initial box size dynamically
        hoomd.run(run_time)
        hoomd.deprecated.dump.xml(group=all, filename="relaxed_" + fileName, all=True)
