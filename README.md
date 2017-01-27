# MorphCT #

The intention of this code is to:

* Parse a coarse-grained morphology (in HOOMD xml format)
* Fine-grain each molecule by returning it to the atomistic representation by building up the polymer chains sequentially from a monomer template
* Output all-atom morphologies in a file that can be read by some sort of structural optimisation/molecular dynamics code to get the rotations of the functional groups correct
* Split it into individual molecules
* Use the fine-grained, correctly-oriented molecules to determine the segments/chromophores that charge carriers are likely to be delocalised along
* Use some kind of quantum chemical method to determine the energy levels of chromophores so that the charge transfer integral can be calculated
* Use the transfer integrals to run a simple Kinetic Monte Carlo algorithm to determine the morphology mobility/charge transport properties

### To run ###

Plop a morphology into inputMorphs (one provided), and call hoomd runMorphCT.py

Note that the majority of the fine-graining process is currently hard-coded to work only with the Huang three-coarse-grained-site P3HT model [PLEASE SEE BRANCH/GENERALIZED TO EXPLORE PROGRESS ON THE GENERALIZED VERSION OF THE CODE DUBBED MORPHCT2.0 THAT CAN BE USED WITH ANY MOLECULES/BLENDS].

### Current Progress ###

* The morphology is being read in, molecules are being fine-grained based on the input template file and rotated according to the locations of the side chains, and a HOOMD xml file being output for the AA morphology. This phase took 3 minutes to complete for a test morphology containing 250 15-mers of P3HT.
* The HOOMD xml file is then being read in and a four-phase HOOMD-Blue simulation is executed:
* 1) A very short simulation for 1E3 timesteps with dt = 1E-5. An nve integrator is used that only permits each atom to move 0.001 distance units per timestep (this is purely to reduce explosions due to atom overlap, hence the tiny runtime). The trajectory is dumped every 10 timesteps, and this phase took 29s to complete for the test morphology.
* 2) A short simulation for 2E3 timesteps with dt = 1E-4. An nvt integrator is used and this phase really only pulls the atoms back that have wandered too far from their structrual equilibrium position during the nve phase. The trajectory is dumped every 10 timesteps, and this phase took 1 minute to complete for the test morphology.
* 3) A long simulation for 1E6 timesteps with dt = 1E-3. An nvt integrator is used and this phase flips the thiophene rings around so that the polymer backbone is ordered correctly. The sidechains are also permitted to move, and tend to flex and stretch to an equilibrium position. Despite the long runtime, this phase polls the kinetic energy every 100 timesteps, and ends the simulation when the kinetic energy reaches the minimum kinetic energy for at least 1000 timesteps (this catches the molecules before they begin to coil and the sidechains wrap around which corresponds to a dip in PE but spike in KE). The trajectory is dumped every 100 timesteps, and this phase took 2 minutes to complete for the test morphology.
* 4) A medium-length simulation for 5E4 timesteps with dt = 1E-3. An nvt integrator is used on the thiophene ring atoms only (i.e. the sidechains are locked in place). This allows the chain backbone to find the lowest energy equilibrium position which is required for the DFT calculations. Equilibrium in this phase is attained quickly, and so the simulation is truncated manually at 50,000 timesteps. The trajectory is dumped every 1E4 timesteps, and this phase took 9 minutes to complete for the test morphology.
* The morphology is then being split into multiple chromophores (currently of length single repeat unit)
* ZINDO/S calculations are performed to determine the electronic transfer integrals between pairs of chromophores
* Kinetic Monte Carlo simulations are then performed using Marcus hopping theory to determine the carrier mobility through the input morphology, which can then be plot as a function of varying state-point value (e.g. annealing temperature).