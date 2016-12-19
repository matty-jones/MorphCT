![alt-text](logo.png "MorphCT Logo")
# MorphCT #

The intention of this code is to:

* Parse any coarse-grained morphology input in the HOOMD xml format. Atomistic morphologies should work too, but the infrastructure is currently set up to map single simulation elements to multiple atoms, so there will likely be plenty of redundant code.
* Fine-grain each morphology by reading in an atomistic molecular template (containing either the full molecule or a single repeat unit in the case of polymer chains) and mapping the atoms to the appropriate simulation elements in the input morphology. The output morphology is geometrically relaxed into a realistic conformation.
* Split the morphology into individual molecules (if required for external analysis e.g. polymer characteristics or quantum chemical computation)
* Determine the electronic chromophores that charge carriers (both electrons and holes for acceptor and donor material respectively) are likely to be delocalised along
* Perform fast quantum chemical calculations (semi-empirical ZINDO/S) to determine the frontier molecular orbitals of each chromophore and calculate the electronic transfer integrals between each neighbouring pair (within a given cut-off)
* Use the transfer integrals to run a simple Kinetic Monte Carlo algorithm to determine the charge carrier mobility, which can be used to infer morphology charge transport properties


### Package Contents and Directory Structures ###

Map of the repo here

### General Code Comments ###

Discussion of data structures and initialisations
I.E. PARAMETER FILES AND HOW TO GENERATE

### Getting Started ###

Plop a morphology into inputMorphs (one provided), and call hoomd runMorphCT.py

Note that the majority of the fine-graining process will have to be hard-coded depending on the polymer chain that we are investigating, so for now it is hard-coded to work only with the Huang three-coarse-grained-site P3HT model.

The fine-graining is handled by code/fineGrainer.py and the execution of hoomd is handled by code/runHoomd.py. code/helperFunctions.py are generic functions/xml parsing subroutines that are shared between the worker modules.

runHoomd.py now fully supports restarting - the output directory is examined to determine the most recently completed stage, and the simulation continues by reading in the most recently output xml file. The fourth phase (the long one) also dumps 100 reset files during its run which can be restarted from if required.

### Specific Module Details ###

Describe the process of how each module works here (and a breakdown of how long it takes)

### Analysis Programs Provided ###

Brief description of how to use each thing

### Future Work ###

* Split the final morphology into individual chains that can be exported as .xyz files for the DFT calculations
* Analyse the chains in the morphology to determine whether the volume can be characterised by a subset of chains (to reduce the number of DFT calculations required)
* Obtain energy levels and transfer integrals
* Monte Carlo all up in this shiznit
