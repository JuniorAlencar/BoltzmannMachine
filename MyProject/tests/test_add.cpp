#include <gtest/gtest.h>
#include "../src/add.h"

// Definir um caso de teste
TEST(AddTest1, PositiveNumbers1) {
    EXPECT_EQ(Add(3, 4), 7);
}

// Definir um caso de teste
TEST(AddTest2, PositiveNumbers2) {
    EXPECT_EQ(Add(4, 4), 8);
}

// Ponto de entrada para o executável de testes
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
