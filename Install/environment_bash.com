#!/bin/bash

# . Bash environment variables and paths to be added to a user's ".bash_profile" file.
# . Some of these values may need modifying (e.g. PEPDICE_SCRATCH and PYTHONPATH).

# . The root of the program.
PEPDICE_ROOT=/home/fernando/programs/pepdice/ ; export PEPDICE_ROOT

# . Package paths.
PEPDICE_BABEL=$PEPDICE_ROOT/Babel                     ; export PEPDICE_BABEL           
PEPDICE_CORE=$PEPDICE_ROOT/Core                       ; export PEPDICE_CORE            
PEPDICE_MOLECULE=$PEPDICE_ROOT/MolecularSystem        ; export PEPDICE_MOLECULE       
#PEPDICE_MOLECULESCRIPTS=$PEPDICE_ROOT/pMoleculeScripts-1.9.0 ; export PEPDICE_PMOLECULESCRIPTS 

# . Additional paths.
PEPDICE_PARAMETERS=$PEPDICE_ROOT/Parameters                                   ; export PEPDICE_PARAMETERS
PEPDICE_SCRATCH=$PEPDICE_ROOT/scratch                                         ; export PEPDICE_SCRATCH   
#PEPDICE_STYLE=$PEPDICE_PARAMETERS/ccsStyleSheets/defaultStyle.css ; export PEPDICE_STYLE     

# . The python path.
PYTHONPATH=$PYTHONPATH:$PEPDICE_ROOT/Babel:$PEPDICE_ROOT/Core:$PEPDICE_ROOT/MolecularSystem ; export PYTHONPATH
