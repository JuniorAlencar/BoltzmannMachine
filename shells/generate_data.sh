#!/bin/bash

# Check for correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <N_spins> <M_states> <seed>"
    echo "You must provide exactly 3 arguments:"
    echo "  N_spins  -> Number of spins"
    echo "  M_states -> Number of states to generate"
    echo "  seed     -> Random seed"
    exit 1
fi

N_spins=$1
M_states=$2
seed=$3

lst=("exact" "metropolis" "parallel_tempering")

# Path to executable
exe_path="../bins/gen_data"

# Loop over methods
for i in "${lst[@]}"
do
    echo "Running: $exe_path $N_spins $M_states $i $seed"
    $exe_path "$N_spins" "$M_states" "$i" "$seed"
done
