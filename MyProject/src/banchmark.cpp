#include "banchmark.hpp"

// Function to perform benchmark and return time and memory usage
std::tuple<long long, double> benchmarkFunction(void (*func)()) {
    // Start measuring time
    auto start = std::chrono::high_resolution_clock::now();

    // Call the function to benchmark
    func();

    // End measuring time
    auto end = std::chrono::high_resolution_clock::now();

    // Compute the duration
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    // Get memory usage
    unsigned long long memoryUsed = 0;
    FILE* file = fopen("/proc/self/statm", "r");
    if (file) {
        unsigned long long size, resident, shared, text, lib, data, dt;
        fscanf(file, "%llu %llu %llu %llu %llu %llu %llu", &size, &resident, &shared, &text, &lib, &data, &dt);
        fclose(file);
        memoryUsed = resident * (unsigned long long)sysconf(_SC_PAGESIZE);
    }

    // Convert memory usage to GB
    double memoryUsedGB = static_cast<double>(memoryUsed) / (1024 * 1024 * 1024);

    // Return the time taken and memory used in GB
    return std::make_tuple(duration.count(), memoryUsedGB);
}