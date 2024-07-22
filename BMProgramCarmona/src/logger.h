// logger.h
// logger.h
#ifndef LOGGER_H
#define LOGGER_H

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/sinks/basic_file_sink.h"

#include <memory>
#include <vector>

// Function to get the global logger
inline std::shared_ptr<spdlog::logger> get_logger() {
    static std::shared_ptr<spdlog::logger> logger = [](){
        // Create console sink
        auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();

        // Create file sink
        auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("logs/global_log.log", true);

        // Combine sinks
        std::vector<spdlog::sink_ptr> sinks {console_sink, file_sink};

        // Create and register logger
        auto combined_logger = std::make_shared<spdlog::logger>("bm", begin(sinks), end(sinks));
        spdlog::register_logger(combined_logger);

        return combined_logger;
    }();

    return logger;
}

#endif // LOGGER_H

