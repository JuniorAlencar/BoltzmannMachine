./ProcessingData $1
./BMfinal $1
./SpecificHeat $1
#./Magnetization_T $1
#./Matriz_Jij $1

path_in='../Results/SeparateData/'
path_out='../Results/Comparativo/'

#Magnetização
paste $path_in'mi-exp/mi_exp_'$1'.dat' $path_in'mi-ising/mi_ising_'$1'.dat' | sed 's/\t/ /g' > $path_out'Magnetização/mag_exp_ising_'$1'.dat'

#Correlação Pearson
paste $path_in'Pij-exp/Pij_exp_'$1'.dat' $path_in'Pij-ising/Pij_ising_'$1'.dat' | sed 's/\t/ /g' > $path_out'Correlação/Pij_exp_ising_'$1'.dat'

#Covariancia
paste $path_in'Cij-exp/Cij_exp_'$1'.dat' $path_in'Cij-ising/Cij_ising_'$1'.dat' | sed 's/\t/ /g' > $path_out'Covariancia/Cij_exp_ising_'$1'.dat'

#sisj
paste $path_in'sisj-exp/sisj_exp_'$1'.dat' $path_in'sisj-ising/sisj_ising_'$1'.dat' | sed 's/\t/ /g' > $path_out'sisj/sisj_exp_ising_'$1'.dat'

#Tijk
paste $path_in'Tijk-exp/Tijk_exp_'$1'.dat' $path_in'Tijk-ising/Tijk_ising_'$1'.dat' | sed 's/\t/ /g' > $path_out'Tripleto/Tijk_exp_ising_'$1'.dat'

#sisjsk
paste $path_in'sisjsk-exp/sisjsk_exp_'$1'.dat' $path_in'sisjsk-ising/sisjsk_ising_'$1'.dat' | sed 's/\t/ /g' > $path_out'sisjsk/sisjsk_exp_ising_'$1'.dat'

