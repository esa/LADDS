/**
 * @file VTUWriter.h
 * @author F. Gratl
 * @date 06.12.21
 */

#pragma once

#include <spdlog/async.h>

#include "breakupModel/util/UtilityContainer.h"
#include "ladds/TypeDefinitions.h"
#include "ladds/distributedMemParallelization/DomainDecomposition.h"
#include "ladds/io/ConfigReader.h"

namespace LADDS {

class VTUWriter {
 public:
  /**
   * Creates a new VTKWriter to a specific file.
   * @param config
   * @param iteration
   * @param decomposition
   */
  explicit VTUWriter(ConfigReader &config, size_t iteration, const DomainDecomposition &decomposition);

  virtual ~VTUWriter();

  /**
   * Write a vtu file (xml style) with the current state of the simulation using the VTKWriter from the breakup code.
   */
  void writeVTU(const AutoPas_t &autopas);

  /**
   * Write a pvtu file which references all vtu files from all ranks in time steps
   * @param iteration
   * @param decomposition
   */
  void writePVTU(ConfigReader &config, size_t iteration, const DomainDecomposition &decomposition);

  /**
   * Write a vtk file in legacy table style (no xml) but with binary storage layout.
   * Without compression this is probably the smallest output file possible.
   * @param iteration
   * @param autopas
   */
  void writeLegacyVTKBinary(size_t iteration, const AutoPas_t &autopas);

  /**
   * Generates a filename for the metadata file (pvtu).
   * @param itstatic eration
   * @return
   */
  std::string filenameMetadata(ConfigReader &config, size_t iteration);

  /**
   * Generates a filename for the payload file (vtu) of a specific rank.
   * @param iteration
   * @param rank If omitted, the current rank number is static chosen.
   * @return
   */
  std::string filenamePayload(ConfigReader &config,
                              size_t iteration,
                              const DomainDecomposition &decomposition,
                              int rank = -1);

 private:
  std::string _loggerName{};

  /**
   * Prints a property of the points.
   * @tparam Property - the type of the property, if it is an array only size = 3 is supported!!!
   * @tparam Data - the class which contains this property
   * @param name - name of the property, e.g. mass
   * @param property - the property (normally a getter of an satellite)
   * @param dataCollection - the data
   */
  template <class Property, class Data, class Container>
  void printProperty(const std::string &name,
                     const std::function<Property(const Data &data)> &property,
                     const Container &dataCollection) const {
    auto logger = spdlog::get(_loggerName);
    if constexpr (std::is_scalar<Property>::value) {
      std::string typeString = [&]() {
        if constexpr (std::is_same_v<Property, float>) {
          return "Float32";
        } else if constexpr (std::is_same_v<Property, double>) {
          return "Float64";
        } else if constexpr (std::is_same_v<Property, int>) {
          return "Int32";
        } else if constexpr (std::is_same_v<Property, long>) {
          return "Int64";
        } else if constexpr (std::is_same_v<Property, unsigned int>) {
          return "UInt32";
        } else if constexpr (std::is_same_v<Property, unsigned long>) {
          return "UInt64";
        } else if constexpr (std::is_same_v<Property, Particle::ActivityState>) {
          return "Int32";
        }
      }();
      logger->info(
          R"(        <DataArray Name="{}" NumberOfComponents="1" format="ascii" type="{}">)", name, typeString);
      for (const auto &date : dataCollection) {
        logger->info("          {}", property(date));
      }
      logger->info(R"(        </DataArray>)");
    } else if constexpr (util::is_stdarray<Property>::value) {
      logger->info(R"(        <DataArray Name="{}" NumberOfComponents="3" format="ascii" type="Float32">)", name);
      for (const auto &date : dataCollection) {
        const auto &array = property(date);
        if (std::size(array) == 3) {
          logger->info("          {} {} {}", array[0], array[1], array[2]);
        }
      }
      logger->info(R"(        </DataArray>)");
    }
  }

  /**
   * Prints the Header of the VTK file.
   * @param size - the number of points
   */
  void printHeader(size_t size) const;

  /**
   * Prints the separator between point properties and point position in cell.
   */
  void printSeparator() const;

  /**
   * Prints the Footer of the VTK file.
   */
  void printFooter() const;
};
}  // namespace LADDS