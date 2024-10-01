#!/bin/bash

# Generate the executables
./compile.sh

# Create folders to Results if don't exist
if [ ! -d "../Results" ] || [ ! -d "../Results_Metropolis" ]; then
	./CreateFolders 
fi

# Run in parallel with arguments in input_multithread.txt
#parallel ./main_single.sh {1} {2} {3} {4} {5} {6} < input_multithread.txt
#parallel ./main_single.sh "$line"
cat input_multithread.txt |parallel --colsep ' ' ./main_single.sh {1} {2} {3} {4} {5} {6}
#cat input_multithread.txt |parallel --colsep ' ' 'echo {1} {2} {3} {4} {5} {6}'
