#!/bin/bash

# GENERATE TO BINS =======================================
# flags to run create folders with libboost
g++ ../scripts/create_folders.cpp -o ../bins/CreateFolders -lboost_filesystem -lboost_system -O3

# Create experimental means
g++ -O3 ../scripts/ProcessingData.cpp -o ../bins/ProcessingData

# Create ising means
g++ -O3 ../scripts/BMfinal.cpp -o ../bins/BMfinal

# Create ising means to specific heat
g++ ../scripts/SpecificHeat.cpp -o ../bins/SpecificHeat

# Create to synthetic_data
g++ ../scripts/gen_data_tests.cpp -o ../bins/gen_data

# Giving permission to execute
chmod 700 ../bins/*

# Create to tests
# g++ ../scripts/tests.cpp -o ../bins/tests

# GENERATE TO FOLDERS RESULTS =====================================

# List of methods
lst=("exact" "metropolis" "parallel_tempering" "swendsen_wang" "wang_landau" "wolff")

# Path to exec
exe_path="../bins/CreateFolders"

# Loop in methods
for i in "${lst[@]}"
do
    echo "Executando $exe_path $i"
    $exe_path "$i"
done
