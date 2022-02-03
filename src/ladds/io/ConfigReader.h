/**
 * @file ConfigReader.h
 * @author F. Gratl
 * @date 20.12.21
 */

#include <autopas/utils/ArrayUtils.h>
#include <autopas/utils/StringUtils.h>

#include <optional>

#include "Logger.h"
#include "yaml-cpp/yaml.h"

#pragma once

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
  ConfigReader(const std::string &configPath, const Logger &logger)
      : config(loadConfig(configPath, logger)), logger(logger){};

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
   *    Has to work with std::to_string() if T != std::string.
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
          // if the value is not a string try to convert it via std::to_string for serialization
          if constexpr (std::is_same_v<T, std::string>) {
            setValue(valuePath, fallbackValue.value());
          } else {
            setValue(valuePath, std::to_string(fallbackValue.value()));
          }
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
   * @param valuePath
   * @param value
   */
  void setValue(const std::string &valuePath, const std::string &value);

 private:
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
   * Map keeping track of all loaded / parsed values.
   * @note This is a sorted map to get a sorted output.
   */
  std::map<std::string, std::string> parsedValues{};
};
