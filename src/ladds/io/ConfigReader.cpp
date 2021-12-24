/**
 * @file ConfigReader.cpp
 * @author F. Gratl
 * @date 20.12.21
 */

#include "ConfigReader.h"

#include <fstream>
#include <iomanip>

YAML::Node ConfigReader::loadConfig(const std::string &cfgFilePath, const Logger &logger) {
  YAML::Node file;
  try {
    file = YAML::LoadFile(cfgFilePath);
  } catch (YAML::Exception &e) {
    logger.log(Logger::Level::warn, "Error loading cfg, loading default file.\n{}", e.msg);
    try {
      file = YAML::LoadFile(defaultCfgPath);
    } catch (YAML::Exception &e) {
      logger.log(Logger::Level::err, "No default file file found. Should be cfg/default_cfg.yaml. Exiting...");
      exit(1);
    }
  }
  return file;
}

void ConfigReader::printParsedValues() {
  int maxKeyLength = 0;
  for (const auto &[key, value] : parsedValues) {
    maxKeyLength = std::max(static_cast<int>(key.size()), maxKeyLength);
  }

  std::stringstream ss{};
  for (const auto &[key, value] : parsedValues) {
    ss << "\n  " << std::left << std::setw(maxKeyLength + 2) << key << ": " << std::right << value;
  }
  SPDLOG_LOGGER_INFO(logger.get(), "Configuration:{}", ss.str());
}

bool ConfigReader::defines(const std::string &valuePath, bool suppressWarning) {
  std::vector<std::string> valuePathVec = autopas::utils::StringUtils::tokenize(valuePath, "/");
  YAML::Node node = Clone(config);
  for (const auto &dir : valuePathVec) {
    // as soon as anything in the path is undefined abort and either use the fallback value or abort.
    if (not node[dir].IsDefined()) {
      if (not suppressWarning) {
        SPDLOG_LOGGER_WARN(logger.get(), "Config value \"{}\" not defined!", valuePath);
      }
      return false;
    }
  }
  return true;
}

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

void ConfigReader::setValue(const std::string &valuePath, const std::string &value) {
  std::vector<std::string> valuePathVec = autopas::utils::StringUtils::tokenize(valuePath, "/");
  setValueAux(config, valuePathVec.begin(), valuePathVec.end(), value);
}

void ConfigReader::dumpConfig(const std::string &filename) const {
  std::ofstream file(filename);
  file << config << std::endl;
  file.close();
}
