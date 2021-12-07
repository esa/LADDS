//
// Created by albert on 05.12.21.
//
#include "Constellation.h"

namespace ConstellationGeneration {

/**
 * generates and returns vector of Constellation objects
 * @param constellationList
 * @param constellationFrequency
 * @param altitudeSpread
 * @return std::pair : Constellation vector , total number of satellites
 */
std::pair<std::vector<Constellation>, size_t> generateConstellations(const std::string &constellationList,
                                                                     size_t constellationFrequency,
                                                                     double altitudeSpread) {
  std::vector<Constellation> constellations;
  size_t satellitesToInsert = 0;
  double altitudeDeviation = altitudeSpread / 3.0;
  if (!constellationList.empty()) {
    auto cons = constellationList;
    int nConstellations = 1;
    for (char con : cons) {
      if (con == ';') {
        nConstellations++;
      }
    }
    for (int i = 0; i < nConstellations; i++) {
      unsigned long offset = cons.find(';', 0);
      if (offset == 0) {
        constellations.emplace_back(Constellation(cons, constellationFrequency, altitudeDeviation));
        break;
      } else {
        constellations.emplace_back(Constellation(cons.substr(0, offset), constellationFrequency, altitudeDeviation));
        cons.erase(0, offset + 1);
      }
    }

    for (int i = 0; i < nConstellations; i++) {
      satellitesToInsert += constellations[i].getConstellationSize();
    }
  }
  return {constellations, satellitesToInsert};
}

}  // namespace ConstellationGeneration
