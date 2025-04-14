# BoltzmannMachine

The main code is present in the BM_program folder. When executing, pay attention only to the shells.

Place the data inside the TidyData folder in comma-separated format, with values ​​1.0/-1.0 without header. Where the columns represent the spins and the lines the states. With **.dat** extension

compile.sh: will generate all executables and all folders to store the results. 
**Run -> ./compile.sh <no arguments>**

generate_data.sh: if you want to generate synthetic data to test the models, run generate_data.sh

**Run -> ./generate_data <N_spins> <M_states> <seed>**
N_spins: Number of spins
M_states: Number of states/samples
seed: random variable

BoltzmannMachine.sh: main code to determine the sets of values ​​J_ij and h_i through the Boltzmann Machine

**Run -> ./BoltzmannMachine.sh <filename> <j_minimum> <h_minimum> <multiply_teq> <multiply_relx> <method>**
**filename**: name of the file from which you want to run the code (without the extension)
**j_minimum**: value decided for the convergence of the values ​​of J_ij
**h_minimum**: value decided for the convergence of the values ​​of h_i
**multiply_teq**: Monte Carlo parameter, where **teq = multiply_teq * N_spins**, where **teq** is the equilibration time.
**multiply_relx**: Monte Carlo parameter, where **relx = multiply_relx * N_spins**, where **relx** is the relaxation time.
**method**: algorithm used, which can be -> exact sum, metropolis or swendsen_wang
