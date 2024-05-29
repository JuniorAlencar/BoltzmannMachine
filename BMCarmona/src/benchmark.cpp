#include "benchmark.h"

#include <fstream>
#include <vector>
#include <string>
#include <numeric>
#include <nlohmann/json.hpp>
#include <experimental/filesystem> // filesystem to c++14
#include <iomanip>  // For std::fixed and std::setprecision
#include <sstream>  // For std::stringstream
#include <iostream> // For std::cerr
#include <unistd.h>
#include <chrono>   // For std::chrono (time process)
#include <functional> // For std::function
#include <tuple> // For std::tuple

using json = nlohmann::json;
namespace fs = std::experimental::filesystem;
using namespace std;

template<typename Func, typename... Args>
bench_values benchmarkFunction(Func&& func, Args&&... args) {
    struct bench_values banch;
    
    // Start measuring time
    auto start = std::chrono::high_resolution_clock::now();

    // Call the function to benchmark
    CallH<Func, Args...>::call(std::forward<Func>(func), std::forward<Args>(args)...);

    // End measuring time
    auto end = std::chrono::high_resolution_clock::now();

    // Compute the duration
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    // Convert duration to milliseconds with 5 decimal places precision
    double timeTaken = static_cast<double>(duration.count()) / 1000.0;

    // Get memory usage
    unsigned long long memoryUsed = 0;
    FILE* file = fopen("/proc/self/statm", "r");
    
    if (file) {
        unsigned long long size, resident, shared, text, lib, data, dt;
        if (fscanf(file, "%llu %llu %llu %llu %llu %llu %llu", &size, &resident, &shared, &text, &lib, &data, &dt) != 7) {
            std::cerr << "Error: Unable to read memory statistics from file." << std::endl;
        }
        fclose(file);
        memoryUsed = resident * (unsigned long long)sysconf(_SC_PAGESIZE);
    }

    // Convert memory usage to GB
    double memoryUsedGB = static_cast<double>(memoryUsed) / (1024 * 1024 * 1024);

    banch.memory = memoryUsedGB;
    banch.time = timeTaken;
    
    // Return the time taken and memory used
    return banch;
}

// Save benchmark values in .json
template<typename Func, typename... Args>
void write_json_benchmark(Func&& func, Args&&... args) {
    // Create folder_benchmark---
    
    // Get the name of the .cpp file
    std::string executableName = fs::path(__FILE__).stem().string();

    // Set the benchmark folder path
    std::string benchmarkFolder = "../benchmark";

    // Create the benchmark folder
    fs::create_directories(benchmarkFolder);
    
    // function_name
    std::string functionName = FUNCTION_NAME(func);

    // Set the benchmark file path
    std::string filename = benchmarkFolder + "/" + executableName + "_" + functionName + ".json";    
    
    // Allocate the memory_usage (Gb) and time_process (ms) vectors
    std::vector<double> memory_usage;
    std::vector<double> time_process;
    
    // Check if the file already exists
    if (fs::exists(filename)) {
        ifstream inputFile(filename);
        if (inputFile.is_open()) {
            json existingData;
            inputFile >> existingData;
            inputFile.close();
            
            if (existingData.contains("memory_usage(Gb)")) {
                memory_usage = existingData["memory_usage(Gb)"].get<std::vector<double>>();
            }
            if (existingData.contains("time_process(ms)")) {
                time_process = existingData["time_process(ms)"].get<std::vector<double>>();
            }
        }
    }

    
    using BenchmarkType = benchmark<int>;  // Replace 'int' with the appropriate type
    // Call the benchmark function with myFunction as argument
    auto result = BenchmarkType::benchmarkFunction(std::forward<Func>(func), std::forward<Args>(args)...);
    
    // Extract variables
    double timeTaken = std::get<0>(result);
    double memoryUsed = std::get<1>(result);
    
    time_process.push_back(timeTaken);
    memory_usage.push_back(memoryUsed);
    
    // Number of implementations
    int n = time_process.size();
    // Sum terms of memory and time
    double sum_time = std::accumulate(time_process.begin(), time_process.end(), 0.0);
    double sum_mem = std::accumulate(memory_usage.begin(), memory_usage.end(), 0.0);
    
    // Avarage of implementations
    double avg_time = sum_time / n;
    double avg_mem = sum_mem / n;

    // Format mean values to 5 decimal places using stringstream
    std::stringstream time_mean_ss;
    std::stringstream memory_mean_ss;
    time_mean_ss << std::fixed << std::setprecision(5) << avg_time;
    memory_mean_ss << std::fixed << std::setprecision(8) << avg_mem;

    // Create JSON object in the desired order
    json j = json::object();
    j["file_test"] = executableName;
    j["method"] = functionName;
    j["memory_mean(Gb)"] = memory_mean_ss.str();
    j["time_mean(ms)"] = time_mean_ss.str();
    j["num_runs"] = n;
    j["memory_usage(Gb)"] = memory_usage;
    j["time_process(ms)"] = time_process;
    
    // Write to the file
    std::ofstream file(filename);
    if (file.is_open()) {
        file << j.dump(4);  // The '4' here is for pretty-printing with an indent of 4 spaces
        file.close();
    } else {
        // Log error or handle the error case
        std::cerr << "Unable to open file for writing: " << filename << std::endl;
        // Handle error, e.g., throw exception or log to a logger
    }
}
