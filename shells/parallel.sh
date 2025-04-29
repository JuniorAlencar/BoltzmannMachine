#!/bin/bash

# Run in parallel with arguments in input_multithread.txt, cat read the file
cat input_multithread.txt |parallel -j 10 --colsep ' ' ./BoltmannMachine.sh {1} {2} {3} {4} {5} {6} {7}