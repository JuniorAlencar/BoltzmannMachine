# $1 -> filename
# $2, $3 -> j_min_erro, h_min_erro (cientific notation, ex: 1.0e-8)
# $4, $5 -> t_eq = n*multi_teq, relx = n*multi_relx
# $6 -> bolean_variable (test), if test=true, create test_files
filename=$1
j_min=$2
h_min=$3
multi_teq=$4
multi_relx=$5
use_exact=$6

./ProcessingData $filename $j_min $h_min $multi_teq $multi_relx $use_exact
./BMfinal $filename $j_min $h_min $multi_teq $multi_relx $use_exact
./SpecificHeat $filename $j_min $h_min $multi_teq $multi_relx $use_exact

#./Magnetization_T $1
#./Matriz_Jij $1

if [ "$use_exact" = "true" ]; then
    path_in='../Results/SeparateData/'
    path_out='../Results/Comparative/'
else
    path_in='../Results_Metropolis/SeparateData/'
    path_out='../Results_Metropolis/Comparative/'
fi

#Magnetização
paste $path_in'mi-exp/mi_exp_'$filename'_err_j_'$j_min'_err_h_'$h_min'_mteq_'$multi_teq'_mrelx_'$multi_relx'.dat' $path_in'mi-ising/mi_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'_mteq_'$multi_teq'_mrelx_'$multi_relx'.dat' | sed 's/\t/ /g' > $path_out'magnetization/mag_exp_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'_mteq_'$multi_teq'_mrelx_'$multi_relx'.dat'

#Correlação Pearson
paste $path_in'Pij-exp/Pij_exp_'$filename'_err_j_'$j_min'_err_h_'$h_min'_mteq_'$multi_teq'_mrelx_'$multi_relx'.dat' $path_in'Pij-ising/Pij_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'_mteq_'$multi_teq'_mrelx_'$multi_relx'.dat' | sed 's/\t/ /g' > $path_out'correlation/Pij_exp_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'_mteq_'$multi_teq'_mrelx_'$multi_relx'.dat'

#Covariancia
paste $path_in'Cij-exp/Cij_exp_'$filename'_err_j_'$j_min'_err_h_'$h_min'_mteq_'$multi_teq'_mrelx_'$multi_relx'.dat' $path_in'Cij-ising/Cij_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'_mteq_'$multi_teq'_mrelx_'$multi_relx'.dat' | sed 's/\t/ /g' > $path_out'covariance/Cij_exp_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'_mteq_'$multi_teq'_mrelx_'$multi_relx'.dat'

#sisj
paste $path_in'sisj-exp/sisj_exp_'$filename'_err_j_'$j_min'_err_h_'$h_min'_mteq_'$multi_teq'_mrelx_'$multi_relx'.dat' $path_in'sisj-ising/sisj_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'_mteq_'$multi_teq'_mrelx_'$multi_relx'.dat' | sed 's/\t/ /g' > $path_out'sisj/sisj_exp_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'_mteq_'$multi_teq'_mrelx_'$multi_relx'.dat'

#Tijk
paste $path_in'Tijk-exp/Tijk_exp_'$filename'_err_j_'$j_min'_err_h_'$h_min'_mteq_'$multi_teq'_mrelx_'$multi_relx'.dat' $path_in'Tijk-ising/Tijk_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'_mteq_'$multi_teq'_mrelx_'$multi_relx'.dat' | sed 's/\t/ /g' > $path_out'triplet/Tijk_exp_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'_mteq_'$multi_teq'_mrelx_'$multi_relx'.dat'

#sisjsk
paste $path_in'sisjsk-exp/sisjsk_exp_'$filename'_err_j_'$j_min'_err_h_'$h_min'_mteq_'$multi_teq'_mrelx_'$multi_relx'.dat' $path_in'sisjsk-ising/sisjsk_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'_mteq_'$multi_teq'_mrelx_'$multi_relx'.dat' | sed 's/\t/ /g' > $path_out'sisjsk/sisjsk_exp_ising_'$filename'_err_j_'$j_min'_err_h_'$h_min'_mteq_'$multi_teq'_mrelx_'$multi_relx'.dat'
