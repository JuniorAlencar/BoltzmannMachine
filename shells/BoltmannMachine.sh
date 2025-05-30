#!/bin/bash

# Check if the correct number of arguments is provided
# File_name j_min_erro h_min_erro

if [ $# -ne 7 ]; then
    echo "Usage: $0 <nome_do_arquivo> <j_min_erro> <h_min_erro> <multi_teq> <multi_relx> <seed> <method>"
    exit 1
fi

filename=$1
j_min=$2
h_min=$3
multi_teq=$4
multi_relx=$5
seed=$6
method=$7

# Check if the input file exists
if [ ! -f "../Data/TidyData/$filename.dat" ]; then
    echo "Error: File '$filename' not found."
    exit 1
fi

# Check if arguments 2 and 3 are valid floats
re_float='^[+-]?([0-9]*[.])?[0-9]+([eE][+-]?[0-9]+)?$'
if ! [[ $j_min =~ $re_float ]]; then
    echo "Error: '$j_min' is not a valid float number."
    exit 1
fi

if ! [[ $h_min =~ $re_float ]]; then
    echo "Error: '$h_min' is not a valid float number."
    exit 1
fi

echo "Running Boltzmann Machine with '$7' method..."

../bins/BMfinal $filename $j_min $h_min $multi_teq $multi_relx $seed $method

sleep 10

echo "Creating comparative files..."

# create files to comparative properties
path_in='../Results/'$7'/SeparateData/'
path_out='../Results/'$7'/Comparative/'

# Magnetization (first moment, si)
paste $path_in'mi-exp/mi_exp_'$1'.dat' $path_in'mi-ising/mi_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat' | sed 's/\t/ /g' > $path_out'magnetization/mag_exp_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat'

# Pearson's correlation
paste $path_in'Pij-exp/Pij_exp_'$1'.dat' $path_in'Pij-ising/Pij_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat' | sed 's/\t/ /g' > $path_out'correlation/Pij_exp_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat'

# Covariance
paste $path_in'Cij-exp/Cij_exp_'$1'.dat' $path_in'Cij-ising/Cij_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat' | sed 's/\t/ /g' > $path_out'covariance/Cij_exp_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat'

# Second moment (sisj)
paste $path_in'sisj-exp/sisj_exp_'$1'.dat' $path_in'sisj-ising/sisj_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat' | sed 's/\t/ /g' > $path_out'sisj/sisj_exp_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat'

# Triplet
paste $path_in'Tijk-exp/Tijk_exp_'$1'.dat' $path_in'Tijk-ising/Tijk_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat' | sed 's/\t/ /g' > $path_out'triplet/Tijk_exp_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat'

# Third momento (sisjsk)
paste $path_in'sisjsk-exp/sisjsk_exp_'$1'.dat' $path_in'sisjsk-ising/sisjsk_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat' | sed 's/\t/ /g' > $path_out'sisjsk/sisjsk_exp_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat'

# hi (random-field)
paste $path_in'hi/hi_real_'$1'.dat' $path_in'hi/hi_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat' | sed 's/\t/ /g' > $path_out'hi/hi_exp_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat'

# Jij (Spin glasses)
paste $path_in'Jij/Jij_real_'$1'.dat' $path_in'Jij/Jij_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat' | sed 's/\t/ /g' > $path_out'Jij/Jij_exp_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat'