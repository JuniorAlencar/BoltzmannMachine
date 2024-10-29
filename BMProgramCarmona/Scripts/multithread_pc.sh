#!/bin/bash

# Define uma função que contêm o código para rodar em paralelo
run_code() {
	time ../bin/bmc ../inputs/$1
}
export -f run_code

arguments=('input_j_min_9.85e-06_h_min_5.81e-05_t_eq_3000_exact.json' 'input_j_min_1.00e-04_h_min_1.00e-04_t_eq_3000_metropolis.json')
parallel run_code :::	 "${arguments[@]}"  
	