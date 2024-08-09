name="sampleN20"

# ----------------------------------------------------------------------------------------------------------------

g++ -O3 ProcessingData.cpp -o ProcessingData3
g++ -O3 BMfinal.cpp -o BMfinal3
g++ -O3 SpecificHeat.cpp -o SpecificHeat3
#g++ -O3 Magnetization_T.cpp -o Magnetization_T
#g++ -O3 Matriz_Jij.cpp -o Matriz_Jij

echo "Processando arquivo "$name" ..."

bash ./Part1_3.sh $name &

cP=$( ps aux | grep "./ProcessingData3" | wc -l )
cB=$( ps aux | grep "./BMfinal3" | wc -l )
cS=$( ps aux | grep "./SpecificHeat3" | wc -l )
#cM=$( ps aux | grep "./Magnetization_T" | wc -l )
#cT=$( ps aux | grep "./Matriz_Jij" | wc -l )

#c=$((cP+cB+cS+cM+CT))
c=$((cP+cB+cS))

while [ $c -ge 17 ] ; do

	cP=$( ps aux | grep "./ProcessingData3" | wc -l )
	cB=$( ps aux | grep "./BMfinal3" | wc -l )
	cS=$( ps aux | grep "./SpecificHeat3" | wc -l )
	#cM=$( ps aux | grep "./Magnetization_T" | wc -l )
	#cT=$( ps aux | grep "./Matriz_Jij" | wc -l )

	#c=$((cP+cB+cS+cM+CT))
	c=$((cP+cB+cS))

	sleep 5
done





