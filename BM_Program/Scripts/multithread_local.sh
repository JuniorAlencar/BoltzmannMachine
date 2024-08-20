#!/bin/bash

	# Define uma função que contêm o código para rodar em paralelo
	run_code() {
		local 	 filename="$1" 
		local 	 err_min_j="$2" 
		local 	 err_min_h="$3" 
		local 	 exact_solution="$4" 
		time ./main_single.sh $filename $err_min_j $err_min_h $exact_solution
		echo "running: $filename $err_min_j $err_min_h $exact_solutions"
 		sleep 2
		}
# Exportar a função usando o módulo Parallel
export -f run_code

args=(('sampleN20' '9.17e-04' '5.76e-04' 'true') ('sampleN20' '9.17e-04' '5.76e-04' 'true') ('sampleN20' '9.17e-04' '5.76e-04' 'true'))

# Executa os comandos em paralelo
for group in "${args[@]}"; do
	clean_group=$(echo "$group" | tr -d '()' | tr -d "'" | xargs)
	echo "$clean_group"
done | parallel -j $(nproc) --colsep ' ' run_code {1} {2} {3} {4} 