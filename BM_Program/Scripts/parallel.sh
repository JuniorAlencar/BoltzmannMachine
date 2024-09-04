#!/bin/bash

# Generate the executables
./compile.sh

# Create folders to Results if don't exist
if [ ! -d "../Results" ] || [ ! -d "../Results_Metropolis" ]; then
	./CreateFolders
fi

# Run in parallel with arguments in input_multithread.txt
parallel --colsep ' ' ./main_single.sh {1} {2} {3} {4} {5} {6} < input_multithread.txt
