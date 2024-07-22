#!/bin/bash

# for i in $( ls *.dat ); do echo $i; done | sed 's/world//g' | awk '{ print "mv world"$1,"World"$1 }' | bash                

# for i in $( ls *.dat ); do echo $i; done | sed 's/configuration//g' | awk '{ print "mv configuration"$1,"World"$1 }' | bash                
# Lista de arquivos
# list_name_file=('configuration2012_N=20' 'configuration2012_N=40' 'configuration2012_N=60' 'configuration2012_N=80'
#                 'configuration2014_N=20' 'configuration2014_N=40' 'configuration2014_N=60' 'configuration2014_N=80'
#                 'configuration2016_N=20' 'configuration2016_N=40' 'configuration2016_N=60' 'configuration2016_N=80'
#                 'configuration2018_N=20' 'configuration2018_N=40' 'configuration2018_N=60' 'configuration2018_N=80')

# list_name_file=('configuration2012_N=20' 'configuration2013_N=20' 'configuration2014_N=20' 'configuration2015_N=20' 'configuration2016_N=20' 'configuration2017_N=20'
#                 'configuration2012_N=40' 'configuration2013_N=40' 'configuration2014_N=40' 'configuration2015_N=40' 'configuration2016_N=40' 'configuration2017_N=40'
#                 'configuration2012_N=60' 'configuration2013_N=60' 'configuration2014_N=60' 'configuration2015_N=60' 'configuration2016_N=60' 'configuration2017_N=60'
#                 'configuration2012_N=80' 'configuration2013_N=80' 'configuration2014_N=80' 'configuration2015_N=80' 'configuration2016_N=80' 'configuration2017_N=80'
#                 'configuration2012_N=100' 'configuration2013_N=100' 'configuration2014_N=100' 'configuration2015_N=100' 'configuration2016_N=100' 'configuration2017_N=100'
#                 'configuration2012_N=120' 'configuration2013_N=120' 'configuration2014_N=120' 'configuration2015_N=120' 'configuration2016_N=120' 'configuration2017_N=120')

# list_name_file=('configuration2012_N=60' 'configuration2013_N=60' 'configuration2014_N=60' 'configuration2015_N=60' 'configuration2016_N=60' 'configuration2017_N=60'
#                 'configuration2012_N=80' 'configuration2013_N=80' 'configuration2014_N=80' 'configuration2015_N=80' 'configuration2016_N=80' 'configuration2017_N=80')

# list_name_file=('World2008_N=60' 'World2008_N=80' 'World2010_N=60' 'World2010_N=80' 'World2012_N=60' 'World2012_N=80'
#                 'configuration2012_N=100' 'configuration2012_N=120' 'configuration2012_N=140' 'configuration2012_N=160'
#                 'configuration2014_N=100' 'configuration2014_N=120' 'configuration2014_N=140' 'configuration2014_N=160'
#                 'configuration2016_N=100' 'configuration2016_N=120' 'configuration2016_N=140' 'configuration2016_N=160'
#                 'configuration2018_N=100' 'configuration2018_N=120' 'configuration2018_N=140' 'configuration2018_N=160')

#list_name_file=$(ls ../Data/TidyData/configuration2014_N=40_sample=*)

# Wolrd samples ------------------------------------------------------------------------------------------------

list_name_pop=$(ls ../Data/TidyData/sampleN20.dat)
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

    bash ./Part1.sh $name  &

    cP=$( ps aux | grep "./ProcessingData" | wc -l )
    cB=$( ps aux | grep "./BMfinal" | wc -l )
    cS=$( ps aux | grep "./SpecificHeat" | wc -l )
    #cM=$( ps aux | grep "./Magnetization_T" | wc -l )
    #cT=$( ps aux | grep "./Matriz_Jij" | wc -l )

    #c=$((cP+cB+cS+cM+CT))
    c=$((cP+cB+cS))

	while [ $c -ge 17 ] ; do

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




