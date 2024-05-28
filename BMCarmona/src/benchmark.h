#ifndef BENCHMARK_H
#define BENCHMARK_H

#include "CallableHelper.h"
#include <utility>
#include <iostream>
#include <vector>
#include <functional>

// Macro to converter function name in string
#define FUNCTION_NAME(func) #func

using namespace std;

// Struct to load memory and time to run process
struct bench_values{
    double memory;
    double time;
};


// Class to calculate and save benchmark values
template<typename Func, typename... Args>
class benchmark{
    private:
        vector<double> memory_usage;
        vector<double> time_process;
        string file_test;
        double memory_mean;
        double time_mean;
        string method;
        int num_run;

    public:
        // Write json_benchmark         
        void write_json_benchmark(Func&& func, Args&&... args);
        
        // Calculate the benchmark values
        bench_values benchmarkFunction(Func&& func, Args&&... args);
        
};

#endif //BENCHMARK_H