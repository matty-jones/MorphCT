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
    return morphologyDict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--temperature", type=float, required=True, help="The temperature of the system")
    parser.add_argument("-p", "--perylothiophene", action="store_true", required=False, help="If present, allow periodic connections to a    dd chromophores to stacks, as well as non-periodic connections (this usually reduces the number of stacks in the system). Defaults to Fa    lse.")
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

        for p in CPAtoms:
            p.mass = 1.000

        for p in CNAtoms:
            p.mass = 0.923

        setCoeffs()

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
        hoomd.run(run_time)
        hoomd.deprecated.dump.xml(group=all, filename="posteql_" + fileName, all=True)
