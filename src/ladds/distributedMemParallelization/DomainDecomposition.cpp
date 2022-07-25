/**
 * @file DomainDecomposition.cpp
 * @author F. Gratl
 * @date 24.05.22
 */
#include "DomainDecomposition.h"

#include <autopas/utils/inBox.h>

std::array<double, 3> LADDS::DomainDecomposition::getGlobalBoxMin() const {
  return globalBoxMin;
}
std::array<double, 3> LADDS::DomainDecomposition::getGlobalBoxMax() const {
  return globalBoxMax;
}
std::array<double, 3> LADDS::DomainDecomposition::getLocalBoxMin() const {
  return localBoxMin;
}
std::array<double, 3> LADDS::DomainDecomposition::getLocalBoxMax() const {
  return localBoxMax;
}
bool LADDS::DomainDecomposition::containsGlobal(const std::array<double, 3> &coordinates) const {
  return autopas::utils::inBox(coordinates, globalBoxMin, globalBoxMax);
}
bool LADDS::DomainDecomposition::containsLocal(const std::array<double, 3> &coordinates) const {
  return autopas::utils::inBox(coordinates, localBoxMin, localBoxMax);
}
autopas::AutoPas_MPI_Comm LADDS::DomainDecomposition::getCommunicator() const {
  return communicator;
}