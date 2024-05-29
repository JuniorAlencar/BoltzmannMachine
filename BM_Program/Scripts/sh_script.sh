

# Wolrd samples ------------------------------------------------------------------------------------------------

#list_name_pop=$(ls ../Data/TidyData/data*)
list_name_pop=$(ls ../Data/TidyData/dataTeste_filter.dat)
#list_name_popB=$(ls ../Data/TidyData/importFortalezaN20.dat)
#list_name_popC=$(ls ../Data/TidyData/importMaracanauN20.dat)


# Soma dos arquivos
list_name_file=("${list_name_pop[@]}")
                #"${list_name_popB[@]}"
                #"${list_name_popC[@]}")

# list_name_file=("${list_name_Municipio_file_2014[@]}" 
#                 "${list_name_Municipio_file_2016[@]}"
#                 "${list_name_Municipio_file_2018[@]}")


# ----------------------------------------------------------------------------------------------------------------

g++ -O3 ProcessingData.cpp -o ProcessingData
g++ -O3 BMfinal.cpp -o BMfinal
g++ -O3 SpecificHeat.cpp -o SpecificHeat
#g++ -O3 Magnetization_T.cpp -o Magnetization_T
#g++ -O3 Matriz_Jij.cpp -o Matriz_Jij

for name in ${list_name_file[@]}; do

    #Amostras
    name=${name#"../Data/TidyData/"}
    name=${name%".dat"}

    echo "Processando arquivo "$name" ..."

    bash ./Part1.sh $name &

    cP=$( ps aux | grep "./ProcessingData" | wc -l )
    cB=$( ps aux | grep "./BMfinal" | wc -l )
    cS=$( ps aux | grep "./SpecificHeat" | wc -l )
    #cM=$( ps aux | grep "./Magnetization_T" | wc -l )
    #cT=$( ps aux | grep "./Matriz_Jij" | wc -l )

    #c=$((cP+cB+cS+cM+CT))
    c=$((cP+cB+cS))

	while [ $c -ge 30 ] ; do

        cP=$( ps aux | grep "./ProcessingData" | wc -l )
        cB=$( ps aux | grep "./BMfinal" | wc -l )
        cS=$( ps aux | grep "./SpecificHeat" | wc -l )
        #cM=$( ps aux | grep "./Magnetization_T" | wc -l )
        #cT=$( ps aux | grep "./Matriz_Jij" | wc -l )

        #c=$((cP+cB+cS+cM+CT))
        c=$((cP+cB+cS))

        sleep 5

    done

done




