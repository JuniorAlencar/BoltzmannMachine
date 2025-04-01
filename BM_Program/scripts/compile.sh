#!/bin/bash

# flags to run create folders with libboost
g++ create_folders.cpp -o ../bins/CreateFolders -lboost_filesystem -lboost_system -O3

# Create experimental means
g++ -O3 ProcessingData.cpp -o ../bins/ProcessingData

# Create ising means
g++ -O3 BMfinal.cpp -o ../bins/BMfinal

# Create ising means to specific heat
g++ SpecificHeat.cpp -o ../bins/SpecificHeat

# Create to synthetic_data
g++ synthetic_data.cpp -o ../bins/synt

# Create to synthetic_data
g++ tests.cpp -o ../bins/tests