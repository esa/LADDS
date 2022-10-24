/**
 * @file DecompositionLogger.cpp
 * @author F. Gratl
 * @date 24.10.22
 */

#include "DecompositionLogger.h"

#include <utility>

LADDS::DecompositionLogger::DecompositionLogger(std::string loggerName,
                                                LADDS::ConfigReader &config,
                                                std::string fileExtension)
    : loggerName(std::move(loggerName)),
      maxDigitsIterations(std::to_string(config.getLastIterationNr()).size()),
      fileExtension(std::move(fileExtension)) {}
