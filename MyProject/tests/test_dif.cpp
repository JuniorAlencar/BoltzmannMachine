#include <gtest/gtest.h>
#include "../src/dif.h"

TEST(DifTest1, NegativeNumbers1) {
    EXPECT_EQ(Dif(3, 4), -1);
}

TEST(DifTest2, NegativeNumbers2) {
    EXPECT_EQ(Dif(4, 4), 0);
}

// Ponto de entrada para o executável de testes
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
