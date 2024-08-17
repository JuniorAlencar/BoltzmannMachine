#!/bin/bash

# Define uma função que contêm o código para rodar em paralelo
run_code() {
	time ../bin/bmc ../inputs/$1
}
export -f run_code

arguments=('input_j_min_9.1e-05_h_min_3.8e-05_t_eq_3000_metropolis.json' 'input_j_min_8.7e-05_h_min_4e-05_t_eq_3000_metropolis.json' 'input_j_min_8.7e-05_h_min_3.8e-05_t_eq_3000_metropolis.json' 'input_j_min_8.7e-05_h_min_4.2e-05_t_eq_3000_metropolis.json' 'input_j_min_8.9e-05_h_min_4e-05_t_eq_3000_metropolis.json' 'input_j_min_8.9e-05_h_min_4.2e-05_t_eq_3000_metropolis.json' 'input_j_min_9.1e-05_h_min_4e-05_t_eq_3000_metropolis.json' 'input_j_min_8.9e-05_h_min_3.8e-05_t_eq_3000_metropolis.json' 'input_j_min_9.1e-05_h_min_4.2e-05_t_eq_3000_metropolis.json')
parallel run_code :::	 "${arguments[@]}"  
	