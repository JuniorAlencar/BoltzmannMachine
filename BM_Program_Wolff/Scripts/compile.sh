#!/bin/bash

# flags to run create folders with libboost
g++ create_folders.cpp -o CreateFolders -lboost_filesystem -lboost_system -O3
# Create experimental means
g++ -O3 ProcessingData.cpp -o ProcessingData
# Create ising means
g++ -O3 BMfinal.cpp -o BMfinal
# Create ising means to specific heat
g++ SpecificHeat.cpp -o SpecificHeat