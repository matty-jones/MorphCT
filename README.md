![alt-text](logo.png "MorphCT Logo" =500x)

# PLEASE NOTE, README IS A WORK IN PROGRESS. #

The intention of this code is to:

* Parse any coarse-grained morphology input in the HOOMD xml format. Atomistic morphologies should work too, but the infrastructure is currently set up to map single simulation elements to multiple atoms, so there will likely be plenty of redundant code.
* Fine-grain each morphology by reading in an atomistic molecular template (containing either the full molecule or a single repeat unit in the case of polymer chains) and mapping the atoms to the appropriate simulation elements in the input morphology. The output morphology is geometrically relaxed into a realistic conformation.
* Split the morphology into individual molecules (if required for external analysis e.g. polymer characteristics or quantum chemical computation)
* Determine the electronic chromophores that charge carriers (both electrons and holes for acceptor and donor material respectively) are likely to be delocalised along
* Perform fast quantum chemical calculations (semi-empirical ZINDO/S) to determine the frontier molecular orbitals of each chromophore and calculate the electronic transfer integrals between each neighbouring pair (within a given cut-off)
* Use the transfer integrals to run a simple Kinetic Monte Carlo algorithm to determine the charge carrier mobility, which can be used to infer morphology charge transport properties

---

# Table of Contents #

* [Package and Directory Structures](#markdown-header-package-and-directory-structures)
* [General Code Comments](#markdown-header-general-code-comments)
* [Job Pipeline](#markdown-header-job-pipeline)
* [Example P3HT Simulation](#markdown-header-example-p3ht-simulation)
* [Analysis Programs Included](#markdown-header-analysis-programs-included)
* [Future Work](#markdown-header-future-work)

---

# Package and Directory Structures #

* ## MorphCT ##
  
  *The main `MorphCT` directory contains all of the required subdirectories, as well as the infrastructure for generating new MorphCT jobs. These jobs are executed by running a parameter file (`parXX.py`), which can be created manually by modifying the template, or by using `generateParameterFile.py`, and answering the questions posed. It is recommended that for each coarse-graining scheme, `generateParameterFile.py` is used only once and then manually copied and modified as desired to account for different inputs, forcefields or simulation parameters to save time.
  The parameter file then calls `runMorphCT.py` to execute the requested MorphCT modules with the parameters outlined in `parXX.py`.
  For your convenience, an example SLURM submission script is also included as `submitMorphXXCT.sh`, which is calibrated to run on Boise State University's Kestrel HPC Cluster.*

  * ### analysis ###
  
    *The `analysis` directory contains several useful scripts for the preparation and analysis of input and output files. A complete description of these and more analysis scripts is given in [Analysis Programs Included](#analysis). Additionally, the user is encouraged to create their own analysis scripts for their system, to utilise the MorphCT data however they see fit.*
  
  * ### code ###
  
    *The `code` directory contains all of the important classes and functions for the various modules that can be executed in MorphCT. A summary of the function of each module and the important features available is given in [Job Pipeline](#pipeline). MorphCT's modular design permits additional modules to be added to improve functionality. The `helperFunctions.py` contains the generic functions/xml parsing subroutines that are shared between the worker modules.*
  
  * ### inputCGMorphs ###
    
    *The `inputCGMorphs` directory contains the coarse-grained simulations that the user would like to input into MorphCT, in order to determine the charge-transport characteristics. The `parXX.py` file selectes the morphology to be considered. The pipeline expects the input morphology to be in HOOMD-Blue XML format. While one of the main features of MorphCT is to fine-grain CG morphologies back to their atomistic representation, there is no reason that the pipeline will not function on an already atomistic morphology, however a significant proportion of code is likely to be redundant, and the user will have to spend significantly more time tuning the `parXX.py` parameter file to work with the infrastructure.*
  
  * ### templates ###
    
    *The `templates` directory contains blank versions of all the files required to make MorphCT run. For instance, `parTemplate.py` is a blank `parXX.py` file that `generateParameterFile.py` reads in and modifies in order to generate a new MorphCT job. `template.inp` is the blank ZINDO/S input file for the ORCA simulation suite, which is used to perform the electronic structure calculations. `sample.sh` is a SLURM script which can be used to send MorphCT jobs to a cluster. Finally, `template.xml` is a blank HOOMD XML file, which is populated with the atomistic morphology after the relaxed atomic positions have been determined.
    Within this directory, the user should place their atomistic templates corresponding to the coarse-grained morphologies in the inputCGMorphs directory. An example has been provided in the repository and explained in [Example - P3HT Simulation](#example).*

---

# General Code Comments #

MorphCT is a modular 'one-size-fits-all' virtual organic electronic device simulator, that crosses all of the important lengthscales for device operation, from the angstroem lengthscale for individual molecules, right up to the hundreds of nanometers required to simulate the organic active layer. In order to effectively transfer information between the standalone modules, a variety of file types and data structures are used:

* .xml (HOOMD)
    * The real workhorse of MorphCT - the input coarse-grained morphologies are entered into the pipeline as HOOMD .xml files, the coarse-grained to fine-grained templates are HOOMD .xml files, and the fine-graining phases produces HOOMD .xml files for visualization purposes.

Discussion of data structures and initialisations
I.E. PARAMETER FILES AND HOW TO GENERATE, PICKLES etc.

Also, don't forget to talk about restarting and how to do it.

---

# Job Pipeline #

* ### Prerequisites ###

1. MorphCT is written in the `python` 2.7 programming language and requires the following modules from the `scipy` stack in addition to the standard python library modules: `numpy`, `matplotlib`, `scipy.sparse`, and `mpl_toolkits` (optional for 3D plotting). Support for Python 3.X will soon be added.
2. In order to use the fine-graining modules `runHoomd` and `extractMol`, the Molecular Dynamics simulation package HOOMD-Blue is required. The source code for HOOMD-Blue can be found at [https://bitbucket.org/glotzer/hoomd-blue/](https://bitbucket.org/glotzer/hoomd-blue/). Note that the current version of MorphCT uses HOOMD 1.3 for its fine-graining procedure and so HOOMD 2.X is currently unsupported. The installation of HOOMD-Blue should automatically set the correct `PATH` and `PYTHONPATH` environment variables, and this can be tested by ensuring that typing `hoomd` into the terminal brings up a python IDE.
3. In order to use the quantum-chemical module `executeZINDO`, the quantum-chemical simulation package ORCA is required. The pre-compiled binaries for ORCA can be downloaded from [https://orcaforum.cec.mpg.de/downloads.php](https://orcaforum.cec.mpg.de/downloads.php). MorphCT has been set up to use ORCA 3.0.3. Previous versions of ORCA may be compatible but are unsupported. After installing ORCA, the `ORCA_BIN` environment variable must be set as the directory that contains the `orca` binary, using the `export` command on Linux-based machines.

* ### Input Files ###

*List input files required and where to put them.*

* ### Generating Parameter Files ###

*Explain everything in the parameter file and how to generate*

* ### Fine-Graining Morphology ###

*For each module be sure to say what the parameter flag is, and how (basically) it works.*

* ### Extracting Molecules ###

* ### Identifying Chromophores ###

* ### Performing Quantum Chemical Calculations ###

* ### Calculating Electronic Transfer Integrals ###

* ### Executing Kinetic Monte Carlo Simulations ###

* ### Analysing Data ###

*Talk about how to extract the mobility data after the KMC*

---

# Example - P3HT Simulation #

For the user's convenience, an example has been provided in the repository. `mid3HT.xml` is an atomistic representation for a single repeat monomer of poly(3-hexylthiophene), which corresponds to the coarse-grained morphologies `p1-L15-f0.0-P0.1-TX.X-e0.5.xml` provided in `inputCGMorphs`. The parameter file `par00.py` shows BLAH BLAH

INCLUDE A BREAKDOWN OF HOW LONG EACH MODULE TAKES

---

# Analysis Programs Included #

For instance, `KMCOut` parses the output KMC pickle file and determines the carrier mobility, anisotropy and connectivity of the simulated morphology. `plotTI` allows for extensive plotting of the electronic properties between statepoints, highlighting the distributions of frontier molecular orbitals, transfer integrals and both intra- and inter-molecular hopping rates. `trimMorphology` can operate on the coarse-grained input file to remove particular coarse-grained elements, while maintaining a consistent intra-molecular constraints list (i.e. remove a CG site type and all bonds, angles, dihedrals and impropers associated with that type while leaving the rest of the morphology intact).    

---

# Future Work #

* Update MorphCT to use Python 3.X
* Update MorphCT to use HOOMD 2.X
* Enhance the Kinetic Monte Carlo simulations to permit full device characterisation (i.e. include exciton generation, multiple carriers, electrical contacts, carrier recombination, separation, and injection, etc.)
* Remove all instances of the `carrierList` in the main morphology pickle as the carriers are treated separately
* Consider removing the `chromophoreList` and treat it separately too?
* Benchmark code and optimise the most commonly-called subroutines

### MJ TO DOs ###

* Merge the generalized branch onto Master when I'm happy everything is working
* Remove redundant analysis scripts and update readme
* Find and fix the mystery seg fault caused by ORCA on Kestrel
* Check which templates we actually use and remove the ones we don't. Then update the readme
