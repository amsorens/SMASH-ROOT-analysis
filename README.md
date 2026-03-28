 .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..
:                                                     :
: Welcome to Agnieszka's SMASH ROOT analysis code!    :
:                                                     :
: Copyright (c) 2026 Agnieszka Sorensen               :
:                                                     :
 .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..

This is an analysis suite for ROOT output files from SMASH simulations.
Currently, the code supports the following analysis types:

SMASH Collider modus:
* flow (basic flow, time evolution of basic flow)
* multiplicity
* yields



::: Directory structure :::

This analysis suite assumes the following structure: a directory "SMASH_results" (can
be any other name) contains a directory called "data" (note: **must** be that name) and
"analysis" (can be any other name). The "data" directory must contain results of SMASH
simulations in directories named with integer numbers, which **must** start from zero
(note that gaps in numbering are allowed, so that {"0", "1", "2", "3", "4", ...} and
{"0", "3", "4", "5", "27", ...} are both fine). The analysis suite must be copied into
"analysis", for example by default as "SMASH-ROOT-analysis" (can be any other name). 
In "SMASH-ROOT-analysis", there is a "src" directory with the source code, the "input"
directory with example config files, and the user creates a "build" directory (can be
any other name) where the analysis suite is compiled. Thus we have
 SMASH_results/                   # Top-level directory (can be any name)
 ├── data/                        # Must be named exactly "data"
 │   ├── 0/                       # Integer-named subdirectories (must start from 0)
 │   ├── 3/
 │   ├── 4/
 │   └── ...
 └── analysis/                    # Can be any name
     └── SMASH-ROOT-analysis/     # Analysis suite (can be any name)
         ├── input/               # Directory with example config file (do NOT remove)
         ├── src/                 # Source code directory (fixed name)
         └── build/               # User-created build directory (can be any name)

Results of analyses can be found in the "SMASH_results/analysis" directory (with some
exceptions; for example, some of the Basic Flow analysis results are also put in a
directory above that, expecting input from different data sets, e.g., for different beam
energies; as another example, some histograms from the Multiplicity analysis are put in
a separate directory within the "analysis" directory to avoid chaos), which is two
directories above the "build" directory.



::: Configuration file :::

Which data directories are taken into account while performing the analysis is defined
in the "Directories" section of the analysis_config.txt file, which after compiling the
code can be found in the "build" directory.
The particular analysis type(s) to be performed are specified in the "Analyses" section
of the config file. The properties of the analysis objects (i.e., pT cuts for flow) are
established either based on information from the config file or take default values from
constants.h and/or analysis object constructors. 



::: Compilation and usage :::

To run the analysis, create an "analysis" directory at the same level as your "data"
directory (see above). Navigate into the "analysis" directory and copy SMASH-ROOT-analysis
there, either directly from GitHub or from wherever you are developing the code[1]. In the
"SMASH-ROOT-analysis" directory, create and navigate into a "build" directory. There,
execute
$ cmake ..
$ make
to compile. Then, specify the range of directories and types of analysis in the
/build/analysis_config.txt" file. To run, execute
$ ./SMASH_ROOT_Analysis

[1] If you develop the SMASH-ROOT-analysis repo and copy the files from there, you will
also copy the hidden git files. These will then ask for a confirmation if you attempt to 
remove the copied directory. To silence these warnings while removing the copy, use the
-rf flag, i.e.,
$ rm -rf SMASH-ROOT-analysis


