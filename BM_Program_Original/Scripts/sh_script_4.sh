name="sampleN20"

# ----------------------------------------------------------------------------------------------------------------

g++ -O3 ProcessingData.cpp -o ProcessingData4
g++ -O3 BMfinal.cpp -o BMfinal4
g++ -O3 SpecificHeat.cpp -o SpecificHeat4
#g++ -O3 Magnetization_T.cpp -o Magnetization_T
#g++ -O3 Matriz_Jij.cpp -o Matriz_Jij

echo "Processando arquivo "$name" ..."

bash ./Part1_4.sh $name &

cP=$( ps aux | grep "./ProcessingData4" | wc -l )
cB=$( ps aux | grep "./BMfinal4" | wc -l )
cS=$( ps aux | grep "./SpecificHeat4" | wc -l )
#cM=$( ps aux | grep "./Magnetization_T" | wc -l )
#cT=$( ps aux | grep "./Matriz_Jij" | wc -l )

#c=$((cP+cB+cS+cM+CT))
c=$((cP+cB+cS))

while [ $c -ge 17 ] ; do

	cP=$( ps aux | grep "./ProcessingData4" | wc -l )
	cB=$( ps aux | grep "./BMfinal4" | wc -l )
	cS=$( ps aux | grep "./SpecificHeat4" | wc -l )
	#cM=$( ps aux | grep "./Magnetization_T" | wc -l )
	#cT=$( ps aux | grep "./Matriz_Jij" | wc -l )

	#c=$((cP+cB+cS+cM+CT))
	c=$((cP+cB+cS))

	sleep 5
done





