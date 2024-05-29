#include <gtest/gtest.h>

#include "../src/exact_solution.h"
#include "../src/experimental_means.h"
#include "../src/forwardmethod.h"
#include "../src/logger.h"
#include "../src/metropolis.h"
#include "../src/network.h"
#include "../src/nr3.h"
#include "../src/read_input_json.h"
#include "../src/read_json.h"
#include "../src/write_json.h"
#include "../src/benchmark.h" // Benchmark calculate

// Ponto de entrada para o executável de testes
// // Definir um caso de teste
// TEST(AddTest1, PositiveNumbers1) {
//     EXPECT_EQ(Add(3, 4), 7);
// }


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}