# ReconstructCDRs

Cloning:
git clone --recursive https://github.com/gabe-raulet/ReconstructCDRs.git

Use the --recursive option because ReconstructCDRs contains a submodule repository for seqan

Compilation:
1. Make a new directory called build.
2. cd into build and type 'cmake .. && make'
3. Use ReconstructCDRs binary executable either by typing full path each time, or set environment PATH variable to include path/to/ReconstructCDRs/build

Usage:
ReconstructCDRs [readFileNameL] [readFileNameR] [bestIsotypesFileName] [organism]

ReconstructCDRs -h for options

organism = 'human' or 'mouse'
