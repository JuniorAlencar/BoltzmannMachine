#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <utility>
using namespace std;

struct bench_values{
    double memory;
    double time;
};

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
    // write json_benchmark 
    template<typename Func, typename Args>
    void write_json_benchmark(Func func, const std::string& functionName, Args args);
    
    // calculate benchmark values
    template<typename Func, typename Args>
    bench_values benchmarkFunction(Func func, const std::string& functionName, Args args);
};

#endif //BENCHMARK_H