/**
 * @file Logger.h
 * @author F. Gratl
 * @date 28.06.21
 */

#pragma once

#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/ostream_sink.h>
#include <spdlog/spdlog.h>

#include <iostream>
#include <string>
#include <utility>

class Logger {
 public:
  explicit Logger(std::string name = "laddsLog", std::ostream &ostream = std::cout);

  ~Logger();

  using Level = spdlog::level::level_enum;

  /**
   * Wrapper around spdlog::log().
   * @tparam Args
   * @param lvl
   * @param fmt
   * @param args
   */
  template <class... Args>
  void log(const Level &lvl, spdlog::string_view_t fmt, const Args &... args) const {
    spdlog::get(_name)->log(lvl, fmt, args...);
  }

  /**
   * Wrapper around spdlog::get().
   * @return Shared pointer to the actual logger.
   */
  [[nodiscard]] auto get() const;

  /**
   * Getter for the logger name.
   * @return Logger name as std::string.
   */
  [[nodiscard]] const std::string &getName() const;

 private:
  std::string _name;
};
