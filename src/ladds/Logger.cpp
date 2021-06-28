/**
 * @file Logger.cpp
 * @author F. Gratl
 * @date 28.06.21
 */

#include "Logger.h"

Logger::Logger(std::string name, std::ostream &ostream) : _name(std::move(name)) {
  auto sink = std::make_shared<spdlog::sinks::ostream_sink_mt>(ostream);
  auto logger = std::make_shared<spdlog::logger>(_name, sink);
  spdlog::register_logger(logger);
}

Logger::~Logger() { spdlog::drop(_name); }

const std::string &Logger::getName() const {
  return _name;
}
