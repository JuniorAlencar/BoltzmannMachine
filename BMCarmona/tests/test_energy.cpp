#include <gtest/gtest.h>
#include "../src/energy.h"    // Energy test
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

