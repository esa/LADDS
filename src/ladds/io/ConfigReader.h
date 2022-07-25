/**
 * @file ConfigReader.h
 * @author F. Gratl
 * @date 20.12.21
 */

#include <autopas/utils/ArrayUtils.h>
#include <autopas/utils/StringUtils.h>
#include <autopas/utils/WrapMPI.h>

#include <optional>

#include "Logger.h"
#include "yaml-cpp/yaml.h"

#pragma once

namespace LADDS {
/**
 * Class that wraps the YAML file, tracks all access into it and can print what was accessed.
 */
class ConfigReader {
 public:
  /**
   * Constructor that loads a YAML file from disk.
   * @param configPath
   * @param logger
   */
  ConfigReader(const std::string &configPath, const Logger &logger, int rank = 0, int numRanks = 1);
  ;

  size_t newParticleID();

  /**
   * Constructor that uses an already loaded YAML tree.
   * @param config
   * @param logger
   */
  ConfigReader(const YAML::Node &config, const Logger &logger) : config(config), logger(logger){};

  /**
   * Primary function to retrieve values from the underlying YAML data.
   * @tparam T Type as which the retrieved data shall be read.
   * @param valuePath Path within the YAML file to the data field.
   * @param fallbackValue Value to use if the field can not be found.
   *    If std::nullopt is used an exception is thrown when no value was found.
   * @param suppressWarning If true no warning is logged when the default value is used.
   * @param fallbackToString If the fallback value is non trivial a function can be passed that provides a string
   * representation.
   * @return Either the value from the YAML file or the fallback value.
   */
  template <class T>
  T get(
      const std::string &valuePath,
      std::optional<T> fallbackValue = std::nullopt,
      bool suppressWarning = false,
      std::function<std::string(T)> fallbackToString = [](const T &t) {
        if constexpr (std::is_same_v<T, std::string>)
          return t;
        else
          return std::to_string(t);
      }) {
    std::vector<std::string> valuePathVec = autopas::utils::StringUtils::tokenize(valuePath, "/");
    YAML::Node node = Clone(config);
    for (const auto &dir : valuePathVec) {
      // as soon as anything in the path is undefined abort and either use the fallback value or abort.
      if (not node[dir].IsDefined()) {
        if (fallbackValue.has_value()) {
          if (not suppressWarning) {
            SPDLOG_LOGGER_WARN(logger.get(),
                               "Config value \"{}\" not defined! Using fallback value: {}",
                               valuePath,
                               fallbackValue.value());
          }
          parsedValues.emplace(valuePath, fallbackToString(fallbackValue.value()));
          // write to value to the datastructure so successive calls will return the same fallback.
          setValue(valuePath, fallbackValue.value());
          return fallbackValue.value();
        } else {
          throw std::runtime_error("Config option not found: " + valuePath);
          return T{};
        }
      }
      node = node[dir];
    }
    parsedValues.emplace(valuePath, node.as<std::string>());
    return node.as<T>();
  }

  /**
   * Checks if a given path is defined in the configuration.
   * @param valuePath
   * @param suppressWarning
   * @return True iff valuePath exists in config.
   */
  bool defines(const std::string &valuePath, bool suppressWarning = false);

  /**
   * Prints all accessed fields and what values were returned.
   */
  void printParsedValues();

  /**
   * Write a YAML file containing the whole configuration object.
   * @param filename
   */
  void dumpConfig(const std::string &filename) const;

  /**
   * Update or insert a value for a given key / value path.
   * If any Node along the path does not exist it will be created.
   * @tparam T Type of the value to write. Has to be std::string, bool, or convertible via std::to_string().
   * @param valuePath
   * @param value
   */
  template <class T>
  void setValue(const std::string &valuePath, const T &value) {
    std::vector<std::string> valuePathVec = autopas::utils::StringUtils::tokenize(valuePath, "/");
    // try to convert value to a string if it isn't already std::string or string literal
    if constexpr (std::is_same_v<T, std::string> or std::is_same_v<T, char[std::extent_v<T>]>) {
      // only set non-empty strings
      if (not static_cast<std::string>(value).empty()) {
        setValueAux(config, valuePathVec.begin(), valuePathVec.end(), value);
      }
    } else if constexpr (std::is_same_v<T, bool>) {
      setValueAux(config, valuePathVec.begin(), valuePathVec.end(), value ? "true" : "false");
    } else {
      setValueAux(config, valuePathVec.begin(), valuePathVec.end(), std::to_string(value));
    }
  }

  /**
   * Getter for the logger.
   * @return logger
   */
  const Logger &getLogger() const;

  /**
   * Determines the number of the fist iteration for the simulation.
   * @return
   */
  size_t getFirstIterationNr();

  /**
   * Determines the number of the last iteration for the simulation.
   * @return
   */
  size_t getLastIterationNr();

 private:
  /**
   * Helper function for setValue recursively iterating through the yaml structure,
   * inserting missing nodes and setting or updating the given value.
   * @tparam T
   * @tparam Iter
   * @param node This is passed by copy as it works similar to a smart pointer.
   * @param begin Start of the array containing the key path.
   * @param end End of the array containing the key path.
   * @param value
   */
  template <typename T, typename Iter>
  void setValueAux(YAML::Node node, Iter begin, Iter end, T value) {
    if (begin == end) {
      return;
    }
    const auto &tag = *begin;
    if (std::next(begin) == end) {
      node[tag] = value;
      return;
    }
    if (!node[tag]) {
      node[tag] = YAML::Node(YAML::NodeType::Map);
    }
    setValueAux(node[tag], std::next(begin), end, value);
  }

  /**
   *  Loads the config file, falls back to default cfg if not found.
   */
  YAML::Node loadConfig(const std::string &cfgFilePath, const Logger &logger);

  /**
   * If no config is given look here for a fallback.
   */
  static constexpr auto defaultCfgPath = "cfg/default_cfg.yaml";

  /**
   * The parsed tree of the YAML file.
   */
  YAML::Node config;

  /**
   * Reference to the logger used for any output.
   */
  const Logger &logger;

  /**
   * Next ID that will be returned the next time newParticleID() is called.
   */
  size_t nextSafeParticleID{};

  /**
   * Last particle ID of the rank-local particle id range which is safe to assign.
   */
  size_t lastSafeParticleID{};

  /**
   * Map keeping track of all loaded / parsed values.
   * @note This is a sorted map to get a sorted output.
   */
  std::map<std::string, std::string> parsedValues{};
};
};  // namespace LADDS