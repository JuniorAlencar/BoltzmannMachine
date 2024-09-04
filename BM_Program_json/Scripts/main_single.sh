#!/bin/bash

# Check if the correct number of arguments is provided
# File_name j_min_erro h_min_erro

echo "Arguments: $1 $2 $3 $4 $5 $6"

if [ $# -ne 6 ]; then
    echo "Usage: $0 <nome_do_arquivo> <j_min_erro> <h_min_erro> <multi_teq> <multi_relx> <use_exact>"
    exit 1
fi

filename=$1
j_min_erro=$2
h_min_erro=$3
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
if ! [[ $j_min_erro =~ $re_float ]]; then
    echo "Error: '$j_min_erro' is not a valid float number."
    exit 1
fi

if ! [[ $h_min_erro =~ $re_float ]]; then
    echo "Error: '$h_min_erro' is not a valid float number."
    exit 1
fi

echo "Processando arquivo "$filename.dat" ..."

# Generate comparative and separate folders with files inside
./Part1_single.sh "$filename" "$j_min_erro" "$h_min_erro" "$multi_teq" "$multi_relx" "$use_exact" &