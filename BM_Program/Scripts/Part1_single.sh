# $1 -> filename
# $2, $3 -> j_min_erro, h_min_erro (cientific notation, ex: 1.0e-8)
# $4 -> bolean_variable (test), if test=true, create test_files

./ProcessingData $1 $4
./BMfinal $1 $2 $3 $4
./SpecificHeat $1 $4
#./Magnetization_T $1
#./Matriz_Jij $1
use_exact=$4
if [ "$use_exact" = "true" ]; then
    path_in='../Results/SeparateData/'
    path_out='../Results/Comparative/'
else
    path_in='../Results_Metropolis/SeparateData/'
    path_out='../Results_Metropolis/Comparative/'
fi

#Magnetização
paste $path_in'mi-exp/mi_exp_'$1'' $path_in'mi-ising/mi_ising_'$1'' | sed 's/\t/ /g' > $path_out'magnetization/mag_exp_ising_'$1''

#Correlação Pearson
paste $path_in'Pij-exp/Pij_exp_'$1'' $path_in'Pij-ising/Pij_ising_'$1'' | sed 's/\t/ /g' > $path_out'correlation/Pij_exp_ising_'$1''

#Covariancia
paste $path_in'Cij-exp/Cij_exp_'$1'' $path_in'Cij-ising/Cij_ising_'$1'' | sed 's/\t/ /g' > $path_out'covariance/Cij_exp_ising_'$1''

#sisj
paste $path_in'sisj-exp/sisj_exp_'$1'' $path_in'sisj-ising/sisj_ising_'$1'' | sed 's/\t/ /g' > $path_out'sisj/sisj_exp_ising_'$1''

#Tijk
paste $path_in'Tijk-exp/Tijk_exp_'$1'' $path_in'Tijk-ising/Tijk_ising_'$1'' | sed 's/\t/ /g' > $path_out'triplet/Tijk_exp_ising_'$1''

#sisjsk
paste $path_in'sisjsk-exp/sisjsk_exp_'$1'' $path_in'sisjsk-ising/sisjsk_ising_'$1'' | sed 's/\t/ /g' > $path_out'sisjsk/sisjsk_exp_ising_'$1''
