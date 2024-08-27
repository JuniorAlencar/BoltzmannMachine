# $1 -> filename
# $2, $3 -> j_min_erro, h_min_erro (cientific notation, ex: 1.0e-8)
# $4 -> bolean_variable (test), if test=true, create test_files
filename=$1
j_min=$2
h_min=$3
use_exact=$4

./ProcessingData $filename $j_min $h_min $use_exact
./BMfinal $filename $j_min $h_min $use_exact
./SpecificHeat $filename $j_min $h_min $use_exact

#./Magnetization_T $1
#./Matriz_Jij $1

if [ "$use_exact" = "true" ]; then
    path_in='../Results/SeparateData/'
    path_out='../Results/Comparative/'
else
    path_in='../Results_Metropolis/SeparateData/'
    path_out='../Results_Metropolis/Comparative/'
fi
paste '../Results_Metropolis/SeparateData/mi-exp/mi_exp_sampleN20_err_j_9.85e-04_err_h_5.81e-04.dat' '../Results_Metropolis/SeparateData/mi-ising/mi_ising_sampleN20_err_j_9.85e-04_err_h_5.81e-04.dat' | sed 's/\t/ /g' > '../Results_Metropolis/Comparative/magnetization/mag_exp_ising_sampleN20_err_j_9.85e-04_err_h_5.81e-04.dat'

#Magnetização
paste $path_in 'mi-exp/mi_exp_'$filename'_err_j_'$j_min'_err_h_'$h_min'.dat' $path_in'mi-ising/mi_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'.dat' | sed 's/\t/ /g' > $path_out'magnetization/mag_exp_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'.dat'

#Correlação Pearson
paste $path_in'Pij-exp/Pij_exp_'$filename'_err_j_'$j_min'_err_h_'$h_min'.dat' $path_in'Pij-ising/Pij_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'.dat' | sed 's/\t/ /g' > $path_out'correlation/Pij_exp_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'.dat'

#Covariancia
paste $path_in'Cij-exp/Cij_exp_'$filename'_err_j_'$j_min'_err_h_'$h_min'.dat' $path_in'Cij-ising/Cij_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'.dat' | sed 's/\t/ /g' > $path_out'covariance/Cij_exp_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'.dat'

#sisj
paste $path_in'sisj-exp/sisj_exp_'$filename'_err_j_'$j_min'_err_h_'$h_min'.dat' $path_in'sisj-ising/sisj_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'.dat' | sed 's/\t/ /g' > $path_out'sisj/sisj_exp_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'.dat'

#Tijk
paste $path_in'Tijk-exp/Tijk_exp_'$filename'_err_j_'$j_min'_err_h_'$h_min'.dat' $path_in'Tijk-ising/Tijk_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'.dat' | sed 's/\t/ /g' > $path_out'triplet/Tijk_exp_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'.dat'

#sisjsk
paste $path_in'sisjsk-exp/sisjsk_exp_'$filename'_err_j_'$j_min'_err_h_'$h_min'.dat' $path_in'sisjsk-ising/sisjsk_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'.dat' | sed 's/\t/ /g' > $path_out'sisjsk/sisjsk_exp_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'.dat'
