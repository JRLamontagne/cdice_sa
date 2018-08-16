#Sensitivity Analysis and Scenario Discovery for DICE

###Overview###
This software package implements the [Dynamic Integrated Climate-Economic model (DICE)] (https://sites.google.com/site/williamdnordhaus/dice-rice) in C++ (CDICE) for scenario discovery and sensitivity analyses purposes.  This repository consists of a C++ implementation of DICE (CDICE), Python supplementary functions and shell scripts to construct the input/output from CDICE for use in sensitivity analyses, and a an example implementation of logistic regression for scenario discovery.  This version of CDICE is adapted from earlier work by [Garner et al.](https://github.com/scrim-network/cdice_doeclim), and has been adapted for use in sensitivity analyses and to run in parallel in HPC settings using MPI.  This repository also relies on the Python-based sensitivity analysis library [SALib](https://github.com/SALib/SALib).  

CDICE is composed of four (4) pieces of code:

CDICE.cpp - The nuts and bolts of the model. CDICEInit.cpp - Function(s) that initialize the DICE model. CDICEMain.cpp - Main function call. Handles input/output of model run. CDICE.h - CDICE header file. Contains structure definitions and function prototypes

To begin, build the model from the source code. The Makefile makes use of mpic++ and should work on Linux (tested on Centos 6.3) or Mac without much trouble. For Windows, use an appropriate C++ MPI compiler. The CDICE executable will except a space delimited file of parameters.  Each line of the file should be a unique scenario, and each line should include 232 values.  The meaning of each of these values is described by parm_desc.txt and an example input is provided in sample_parms.txt.

Next, download or clone [SALib](https://github.com/SALib/SALib) from the Github and place it in the main project directory. An example submission script (using PBS in an HPC setting) that conducts a sobol sensitivity analysis on various output from CDICE and replicates the experiment is included in the Bootstrap directory: Example_Run.sh.  The script begins by reading a configuration file and policy configuration file.  The configuration file indicates which exogenous DICE variables are to be considered in the analysis and their sampling ranges.  The policy configuration file indicates what policy functional form is to be used and how that form should be parameterized.  Examples of both configuration files are included in the repository's main directory.  Next the configuration files are processed to be in the correct form for SALib's Saltelli sampling function, which in turn generates scenario samples, which are in turn configured in the correct format for CDICE.  Each scenario is evaluated using CDICE run in parallel, and the output is re-formmatted to be compatible with SALib's Sobol sensitivity analysis routines.  The remainder of the Example_Run.sh script computes the Sobol Sensitivity Indicies for various CDICE outputs.

The logitModel directory contains an example script for applying logistic regression for scenario discovery using the statsmodels Python package.

The data directory contains all of the data used to conduct the analysis and generate the figures in our paper.

**Requirements:**
Python 2.7.5 including packages: [Numpy](http://www.numpy.org/), [Scipy](https://www.scipy.org/), [statsmodels](https://www.statsmodels.org/stable/index.html), and [mpi4py](https://pypi.org/project/mpi4py/)

**Citation:**
Lamontagne, J.R., P. M. Reed, G. Marangoni, K. Keller and G. G. Garner. Robust pathways to tolerable climate futures require immediate
global action, Nature Climate Change, (in Review), 2018.