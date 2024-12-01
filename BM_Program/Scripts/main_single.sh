#!/bin/bash

# Check if the correct number of arguments is provided
# File_name j_min_erro h_min_erro

if [ $# -ne 6 ]; then
    echo "Usage: $0 <nome_do_arquivo> <j_min_erro> <h_min_erro> <multi_teq> <multi_relx> <use_exact>"
    exit 1
fi

filename=$1
j_min=$2
h_min=$3
multi_teq=$4
multi_relx=$5
use_exact=$6

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

echo "Processando arquivo "$filename.dat" ..."

# Generate comparative and separate folders with files inside
./CreateFolders &
./ProcessingData $filename $use_exact &
./BMfinal $filename $j_min $h_min $multi_teq $multi_relx $use_exact &

# pause for 5 seconds
sleep 5

# create files to comparative properties

if [ "$use_exact" = "true" ]; then
    path_in='../Results/SeparateData/'
    path_out='../Results/Comparative/'
else
    path_in='../Results_Metropolis/SeparateData/'
    path_out='../Results_Metropolis/Comparative/'
fi

# Magnetization (first moment, si)
paste $path_in'mi-exp/mi_exp_'$1'.dat' $path_in'mi-ising/mi_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat' | sed 's/\t/ /g' > $path_out'Magnetização/mag_exp_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat'

# Pearson's correlation
paste $path_in'Pij-exp/Pij_exp_'$1'.dat' $path_in'Pij-ising/Pij_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat' | sed 's/\t/ /g' > $path_out'Correlação/Pij_exp_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat'

# Covariance
paste $path_in'Cij-exp/Cij_exp_'$1'.dat' $path_in'Cij-ising/Cij_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat' | sed 's/\t/ /g' > $path_out'Covariancia/Cij_exp_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat'

# Second moment (sisj)
paste $path_in'sisj-exp/sisj_exp_'$1'.dat' $path_in'sisj-ising/sisj_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat' | sed 's/\t/ /g' > $path_out'sisj/sisj_exp_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat'

# Triplet
paste $path_in'Tijk-exp/Tijk_exp_'$1'.dat' $path_in'Tijk-ising/Tijk_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat' | sed 's/\t/ /g' > $path_out'Tripleto/Tijk_exp_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat'

# Third momento (sisjsk)
paste $path_in'sisjsk-exp/sisjsk_exp_'$1'.dat' $path_in'sisjsk-ising/sisjsk_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat' | sed 's/\t/ /g' > $path_out'sisjsk/sisjsk_exp_ising_'$1'_err_j_'$2'_err_h_'$3'_mteq_'$4'_mrelx_'$5'.dat'


#./SpecificHeat $filename $j_min $h_min $multi_teq $multi_relx $use_exact &




