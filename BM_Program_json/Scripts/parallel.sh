./compile.sh

parallel ./main_single.sh {1} {2} {3} {4} {5} {6} < input_multithread.txt
#parallel echo {1} {2} {3} {4} {5} {6} :::: input_multithread.txt