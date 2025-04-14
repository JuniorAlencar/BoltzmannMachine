# BoltzmannMachine

# BoltzmannMachine

The main code is located in the `BM_program` folder. Execution is performed through the shell scripts provided.  

Place your data inside the `TidyData` folder in **comma-separated `.dat` files** using values `1.0` and `-1.0`, without a header.  
Each column corresponds to a spin, and each row to a state.

---

## 📦 Executables

The main shell scripts provided are:

- `compile.sh`: Compiles all necessary executables and creates the folder structure to store results.  
- `generate_data.sh`: (Optional) Generates synthetic data for testing.
- `BoltzmannMachine.sh`: Main script to compute `J_ij` and `h_i` using the Boltzmann Machine approach.

---

## 🚀 Running Code

```bash
# Compile all executables and create folders for results
./compile.sh



# Generate synthetic data (optional)
# Usage:
# ./generate_data <N_spins> <M_states> <seed>
# Example:
./generate_data 100 5000 42


# Run the Boltzmann Machine to estimate J_ij and h_i
# Usage:
# ./BoltzmannMachine.sh <filename> <j_minimum> <h_minimum> <multiply_teq> <multiply_relx> <method>
#
# Parameters:
# filename       -> name of the input file (without .dat extension)
# j_minimum      -> minimum error for convergence of J_ij (e.g., 1e-5)
# h_minimum      -> minimum error for convergence of h_i  (e.g., 1e-4)
# multiply_teq   -> equilibration time factor: teq = multiply_teq * N_spins
# multiply_relx  -> relaxation time factor: relx = multiply_relx * N_spins
# method         -> algorithm: exact, metropolis, swendsen_wang, parallel_tempering, wang_landau
#
# Example:
./BoltzmannMachine.sh sample_data 1e-5 1e-4 150 2 metropolis