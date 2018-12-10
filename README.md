[![alt-text](morphCTLogo.png "MorphCT Logo" =500x)](https://bitbucket.org/cmelab/morphct)

[![DOI](https://zenodo.org/badge/100152281.svg)](https://zenodo.org/badge/latestdoi/100152281)
[![Run Status](https://api.shippable.com/projects/5ae1ee243e7d2b0700153be4/badge?branch=dev)](https://app.shippable.com/bitbucket/cmelab/morphct) 
[![Coverage Status](https://codecov.io/bb/cmelab/morphct/branch/dev/graph/badge.svg)](https://codecov.io/bb/cmelab/morphct)

# Intellectual Property #

MorphCT has been released under a GPL3 license (please see LICENSE.TXT for conditions). Please read and cite our [Molecular Simulations paper](https://doi.org/10.1080/08927022.2017.1296958), and our [zenodo DOI](https://zenodo.org/badge/latestdoi/100152281) if you use our code in your published work.

---

# Code Description #

The intention of this code is to form a modular, computational pipeline that can powerfully and coherently relate organic molecular morphology on the Angstroem lengthscale to electronic device performance over hundreds of nanometers.

MorphCT accomplishes this by:

* Converting any coarse-grained organic morphology to the atomistic representation for submission to quantum chemical calculations through a process called fine-graining.
* Splitting the morphology into electronically-active chromophores that charge carriers (holes in the case of donor materials, electrons for acceptors) are likely to be delocalised along, and can perform quantized charge hops between.
* Performing high-throughput, fast quantum chemical calculations (QCCs) to obtain the energetic landscape caused by conformational disorder, as well as electronic transfer integrals between chromophore pairs.
* Using these calculated electronic properties as inputs into a kinetic Monte Carlo (KMC) algorithm to simulate the motion of charge carriers throughout the device, allowing carrier mobilities to be obtained (a good proxy for device performance)
* Using the charge transport properties of previous morphological moieties to simulate the various mechanisms as play in an organic photovoltaic device (applied voltage bias, photoinjection, dark-current injection, exciton transport, exciton dissociation, carrier transport, extraction, geminate recombination, bimolecular recombination) in order to generate J-V curves and calculate OPV device efficiencies to guide manufacturing.

---

# Tagged Releases #

* v3.1: Many quality-of-life improvements, issue resolutions, and bug fixes including but not limited to: the ability to write ORCA jobs to a RAM disk for 30x faster calculations, significantly more unit tests to increase code coverage to 75%, full testing the device simulations, better python formatting, and migration from Coveralls to Codecov.io.
* v3.0: MorphCT codebase brought up to PEP8 standard (with line limit of 120 characters), refactored to work as a package, added extensive unit tests to the pipeline (pytest) and added continuous integration support using Shippable and Coveralls.
* v2.2: Additional funcionality added for blend morphologies, multiple donor/acceptor species, variable-range hopping, Voronoi neighbourhood analysis and other performance tweaks. Results in preparation to be submitted in Q2 2018.
* v2.1: MorphCT updated from python 2.7 to python 3.5
* v2.0: Hardcode removed, utility of MorphCT expanded to any organic molecule with more customizable features (support for small molecules included)
* v1.0: MorphCT released under GPL. Hardcoded mobility results for P3HT, results published in [Molecular Simulation](https://doi.org/10.1080/08927022.2017.1296958) 

---

# Getting Started #

## Installation ##

* MorphCT requires Python 3.
* Additionally, due to the number of non-core python modules used in the codebase, it is recommended to use an environment manager to set up and run MorphCT. Prerequisites are listed below if the user desires manual setup.
    * Support is provided for using the Conda environment manager, and users are recommended to download and install the appropriate distribution of [Miniconda](https://conda.io/miniconda.html) for their operating system.
    * After Miniconda has been installed and added to the `$PATH` variable, the environment can be set up by using `conda env create -f environment.yml` from inside MorphCT's root directory. This will use Conda and Pip to install all of the required modules to run MorphCT successfully.
    * The created environment can be activated using `source activate morphct`. This environment will need to be active in order for MorphCT to run.
    * Now that the dependencies are installed, the MorphCT package must be installed with pip by using `pip install -e .`.
* In order to perform the molecular dynamics simulations, the user will require v1.3 of [HOOMD-Blue](http://glotzerlab.engin.umich.edu/hoomd-blue/) to be installed on their system. This can be obtained from the above link, and is freely available for all users. Please read and cite the following papers [#1](https://doi.org/10.1016/j.jcp.2008.01.047), [#2](https://doi.org/10.1016/j.cpc.2015.02.028) when using HOOMD-Blue.
    * The Miniconda environment discussed above includes a build of HOOMD-Blue v1.3, which suffices for the unit tests and small systems (< 1000 particles). Depending on the infrastructure available on the user's system, it may be beneficial to remove this instantiation of HOOMD, and instead use a custom-compiled build to utilise any graphical processing units for larger systems.
    * Note that MorphCT does not yet support v2.X of HOOMD, due to technical issues concerning the implementation of rigid bodies.
* In order to perform the quantum chemical calculations, the user will require [orca](https://cec.mpg.de/orcadownload/index.php) to be installed on their system. Orca is available from the above link and is freely available for academic use. Please read and cite the following [paper](https://doi.org/10.1002/wcms.81) when using orca.
    * Support is officially provided for orca v3.X, but orca v4.X also functions correctly.
    * Note that MorphCT uses distutils to find the orca executable in the path. In the event that the user has multiple orca binaries, this search can be overriden by setting the `ORCA_BIN` environment variable to the quantum chemical suite's binary.

## Dependencies ##

* Required:
    * Python == 3.5
    * Scipy == 0.19.1
    * Matplotlib == 2.0.2

* Optional:
    * HOOMD-Blue == 1.3
    * Orca >= 3.0
    * Pytest
    * Pytest-Cov
    * PyYAML
    * Coveralls
    * ImageMagick

## Running Tests ##

* 300 unit tests are provided to ensure that MorphCT has been set up correctly.
* These can be executed by running `pytest` in the MorphCT root directory, or in `morphct/tests`.

---

# Running Jobs #

## Inputs ##

* After setup, MorphCT requires two types of inputs: the input molecular system (coarse-grained or already atomistic), and a relevant parameter dictionary describing the modules to execute and the parameters to use.
* The input molecular system consists of several xml files depending on the users requirements:
    * Input morphology xml: The main molecular morphology to be input to the pipeline. This should be in HOOMD xml format and can be coarse-grained (in which case, the fine-graining module must be called first before the transport properties can be obtained), or atomistic (if only the transport properties of the system are required). If providing an atomistic morphology, then this is the only molecular system input required.
    * Atomistic Template: If fine-graining is desired, this HOOMD xml file describes the atomstic template to be mapped onto the input coarse-grained system. Examples are given in the `templates` directory for a variety of organic molecule types: `P3HT.xml`, `C60.xml`, `BDTTPD.xml`, `perylene.xml`, `perylothiophene.xml` and `PCBM.xml`.
    * Atomistic Forcefields: If fine-graining is desired, this xml file contains the forcefield parameters describing the atomistic interactions used to relax the morphology to a realistic conformation after fine-graining. Examples are given in the `templates` directory for several of the example molecules `FFP3HT.xml`, `FFC60.xml`, `FFPerylene.xml`, `FFPerylothiophene.xml`.
* An example is provided for the user in `morphct/templates/par_template.py`. This file is heavily documented to describe the functions and features of MorphCT, and should be the first port-of-call for seeing what is possible with the program.
* For device simulations, a set of molecular systems already run through the pipeline (up to the mobility KMC phase - i.e. morphology fine-grained and relaxed, chromophores detected, and energy levels and transfer integrals calculated) for every component moiety of the input device must be present and have its location identified in the parameter file. An additional device input file identifying the arrangement of the morphology moieties in the device must also be present and have its location identified in the parameter file.

## Job Execution ##

* Jobs can be invoked by using the HOOMD v1.3 python wrapper around the relevant parameter file for the job: `hoomd par_template.py`.
* MorphCT has also been heavily tested on several High Performance Computing Clusters (Kestrel, Fry, R1, R2, Comet, Bridges, XStream), and found to scale well. Child jobs are spawned for parallelisable modules (quantum chemical calculations, mobility KMC and device KMC) in the pipeline, so simply submitting the `hoomd par_template.py` command to the HPC resource management system should be sufficient to parallelise the jobs where possible.
* For large systems (> 1000 chromophores), it is recommended to perform the QCC calculations in memory, or at least on a fast file system with high file i/o bandwidth (Orca reads and writes a large number of temporary files to the disk for every calculation, which can quickly clog up a slower file system), by setting the `orca_output_directory` parameter.

---

# Contributing #

* Please note that the [github repository](https://github.com/matty-jones/MorphCT) for MorphCT is a mirror of the codebase stored on [Bitbucket](https://bitbucket.org/cmelab/morphct). Issues and PRs should be raised on Bitbucket (free account).
* Please feel free to [fork](https://confluence.atlassian.com/bitbucket/forking-a-repository-221449527.html) the repository and submit your [pull requests](https://www.atlassian.com/git/tutorials/making-a-pull-request) to the `dev` branch.
* Alternatively, raising issues and voting for the issues most important to you is highly recommended and appreciated.
* All contributions should be PEP8 compliant (comment lines limited to 80 characters, code lines limited to 120 characters).
* All pull requests require approval from at least one reviewer, a successful build on the latest commit to the fork, and no failed builds on the latest commit to the fork to be accepted.
* Please note that if the pull request of your fork includes changes to the shippable.yml, the CI will still use the yml on the destination branch (usually dev). To resolve this, please activate CI for your fork first to ensure build stability, and link the CI results in your pull request.

---

# Maintainers #

* Matthew Jones (mattyjones@boisestate.edu)
* Mike Henry (mikehenry@boisestate.edu)
