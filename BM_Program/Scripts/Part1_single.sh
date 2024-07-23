# $1 -> filename
# $2, $3 -> j_min_erro, h_min_erro (cientific notation, ex: 1.0e-8)
# $4 -> bolean_variable (test), if test=true, create test_files

./ProcessingData $1
./BMfinal $1 $2 $3 $4
./SpecificHeat $1
#./Magnetization_T $1
#./Matriz_Jij $1

test=$4

path_in='../Results/SeparateData/'
path_out='../Results/Comparativo/'

# If test = True, create test folders with multiplies files
if [ "$test" = "true" ]; then
    #Magnetização
    paste $path_in'mi-exp/tests/mi_exp_'$1'_j_'$2'_h_'$3'.dat' $path_in'mi-ising/tests/mi_ising_'$1'_j_'$2'_h_'$3'.dat' | sed 's/\t/ /g' > $path_out'Magnetização/tests/mag_exp_ising_'$1'_j_'$2'_h_'$3'.dat'

    #Correlação Pearson
    paste $path_in'Pij-exp/tests/Pij_exp_'$1'_j_'$2'_h_'$3'.dat' $path_in'Pij-ising/tests/Pij_ising_'$1'_j_'$2'_h_'$3'.dat' | sed 's/\t/ /g' > $path_out'Correlação/tests/Pij_exp_ising_'$1'_j_'$2'_h_'$3'.dat'

    #Covariancia
    paste $path_in'Cij-exp/tests/Cij_exp_'$1'_j_'$2'_h_'$3'.dat' $path_in'Cij-ising/tests/Cij_ising_'$1'_j_'$2'_h_'$3'.dat' | sed 's/\t/ /g' > $path_out'Covariancia/tests/Cij_exp_ising_'$1'_j_'$2'_h_'$3'.dat'

    #sisj
    paste $path_in'sisj-exp/tests/sisj_exp_'$1'_j_'$2'_h_'$3'.dat' $path_in'sisj-ising/tests/sisj_ising_'$1'_j_'$2'_h_'$3'.dat' | sed 's/\t/ /g' > $path_out'sisj/tests/sisj_exp_ising_'$1'_j_'$2'_h_'$3'.dat'

    #Tijk
    paste $path_in'Tijk-exp/tests/Tijk_exp_'$1'_j_'$2'_h_'$3'.dat' $path_in'Tijk-ising/tests/Tijk_ising_'$1'_j_'$2'_h_'$3'.dat' | sed 's/\t/ /g' > $path_out'Tripleto/tests/Tijk_exp_ising_'$1'_j_'$2'_h_'$3'.dat'

    #sisjsk
    paste $path_in'sisjsk-exp/tests/sisjsk_exp_'$1'_j_'$2'_h_'$3'.dat' $path_in'sisjsk-ising/tests/sisjsk_ising_'$1'_j_'$2'_h_'$3'.dat' | sed 's/\t/ /g' > $path_out'sisjsk/tests/sisjsk_exp_ising_'$1'_j_'$2'_h_'$3'.dat'

# Else: create final file
else
    #Magnetização
    paste $path_in'mi-exp/mi_exp_'$1'' $path_in'mi-ising/mi_ising_'$1'' | sed 's/\t/ /g' > $path_out'Magnetização/mag_exp_ising_'$1''

    #Correlação Pearson
    paste $path_in'Pij-exp/Pij_exp_'$1'' $path_in'Pij-ising/Pij_ising_'$1'' | sed 's/\t/ /g' > $path_out'Correlação/Pij_exp_ising_'$1''

    #Covariancia
    paste $path_in'Cij-exp/Cij_exp_'$1'' $path_in'Cij-ising/Cij_ising_'$1'' | sed 's/\t/ /g' > $path_out'Covariancia/Cij_exp_ising_'$1''

    #sisj
    paste $path_in'sisj-exp/sisj_exp_'$1'' $path_in'sisj-ising/sisj_ising_'$1'' | sed 's/\t/ /g' > $path_out'sisj/sisj_exp_ising_'$1''

    #Tijk
    paste $path_in'Tijk-exp/Tijk_exp_'$1'' $path_in'Tijk-ising/Tijk_ising_'$1'' | sed 's/\t/ /g' > $path_out'Tripleto/Tijk_exp_ising_'$1''

    #sisjsk
    paste $path_in'sisjsk-exp/sisjsk_exp_'$1'' $path_in'sisjsk-ising/sisjsk_ising_'$1'' | sed 's/\t/ /g' > $path_out'sisjsk/sisjsk_exp_ising_'$1''
fi