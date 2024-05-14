#include <gtest/gtest.h>
#include "../src/my_code.h"

// Definir um caso de teste
TEST(AddTest, PositiveNumbers) {
    EXPECT_EQ(Add(3, 4), 7);
}

TEST(AddTest, NegativeNumbers) {
    EXPECT_EQ(Add(-3, -4), -7);
}

// Ponto de entrada para o executável de testes
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
