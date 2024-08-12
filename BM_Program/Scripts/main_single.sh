#!/bin/bash

# flags to run create folders with libboost
g++ ProcessingData.cpp -o ProcessingData -lboost_filesystem -lboost_system -O3 
g++ -O3 BMfinal.cpp -o BMfinal
g++ -O3 SpecificHeat.cpp -o SpecificHeat
#g++ -O3 Magnetization_T.cpp -o Magnetization_T
#g++ -O3 Matriz_Jij.cpp -o Matriz_Jij

# Check if the correct number of arguments is provided

# File_name j_min_erro h_min_erro
if [ $# -ne 4 ]; then
    echo "Usage: $0 <nome_do_arquivo> <j_min_erro> <h_min_erro> <use_exact>"
    exit 1
fi

filename=$1
j_min_erro=$2
h_min_erro=$3
use_exact=$4

# Check if the input file exists
if [ ! -f "../Data/TidyData/$filename" ]; then
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

echo "Processando arquivo "$filename" ..."

# Generate comparative and separate folders with files inside
#bash ./Part1_single.sh "$filename" "$j_min_erro" "$h_min_erro" &
./Part1_single.sh "$filename" "$j_min_erro" "$h_min_erro" "$use_exact"&


count_processes() {
    local process_name=$1
    ps aux | grep "./$process_name" | grep -v "grep" | wc -l
}

# Count the running processes initially
cP=$(count_processes "ProcessingData")
cB=$(count_processes "BMfinal")
cS=$(count_processes "SpecificHeat")
#cM=$( ps aux | grep "./Magnetization_T" | wc -l )
#cT=$( ps aux | grep "./Matriz_Jij" | wc -l )

#c=$((cP+cB+cS+cM+CT))
c=$((cP+cB+cS))

# Loop while the total number of processes is greater than or equal to 17
while [ $c -ge 17 ]; do
    cP=$(count_processes "ProcessingData")
    cB=$(count_processes "BMfinal")
    cS=$(count_processes "SpecificHeat")
    # cM=$(count_processes "Magnetization_T")
    # cT=$(count_processes "Matriz_Jij")

    # Total running processes
    # c=$((cP + cB + cS + cM + cT))
    c=$((cP + cB + cS))

    sleep 5
done

# Run the scripts with the three arguments
# bash ./ProcessingData $filename &
# bash ./BMfinal "$filename" "$j_min_erro" "$h_min_erro" "$test" &
# bash ./SpecificHeat "$filename" &
./ProcessingData $filename "$use_exact" &
./BMfinal "$filename" "$j_min_erro" "$h_min_erro" "$use_exact" &
./SpecificHeat "$filename" "$use_exact" &
# bash ./Magnetization_T "$filename" "$j_min_erro" "$h_min_erro" &
# bash ./Matriz_Jij "$filename" "$j_min_erro" "$h_min_erro" &






