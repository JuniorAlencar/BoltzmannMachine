#ifndef BANCHMARK_H
#define BANCHMARK_H

#include <iostream>
#include <chrono>
#include <vector>
#include <tuple>
#include <unistd.h> // Include for sysconf

std::tuple<long long, double> benchmarkFunction(void (*func)());

#endif // DIF_H