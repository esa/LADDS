/**
 * @file ParticleCommunicator.cpp
 * @author F. Gratl
 * @date 01.06.22
 */

#include "ParticleCommunicator.h"

#include "Serialization.h"

autopas::AutoPas_MPI_Request LADDS::ParticleCommunicator::sendParticles(
    const std::vector<LADDS::Particle>::iterator &particlesBegin,
    const std::vector<LADDS::Particle>::iterator &particlesEnd,
    const int receiver,
    autopas::AutoPas_MPI_Comm const &communicator) {
  std::vector<char> serializedParticles;
  Serialization::serializeParticles(particlesBegin, particlesEnd, serializedParticles);
  autopas::AutoPas_MPI_Request sendRequest{};
  autopas::AutoPas_MPI_Isend(serializedParticles.data(),
                             static_cast<int>(serializedParticles.size()),
                             AUTOPAS_MPI_CHAR,
                             receiver,
                             0,
                             communicator,
                             &sendRequest);

  return sendRequest;
}

std::vector<LADDS::Particle> LADDS::ParticleCommunicator::receiveParticles(
    const int sender, autopas::AutoPas_MPI_Comm const &communicator) {
  std::vector<char> serializedParticles;

  autopas::AutoPas_MPI_Status status;
  autopas::AutoPas_MPI_Probe(sender, 0, communicator, &status);

  int receiveBufferSize = 0;
  autopas::AutoPas_MPI_Get_count(&status, AUTOPAS_MPI_CHAR, &receiveBufferSize);
  serializedParticles.resize(receiveBufferSize);

  autopas::AutoPas_MPI_Recv(serializedParticles.data(),
                            receiveBufferSize,
                            AUTOPAS_MPI_CHAR,
                            sender,
                            0,
                            communicator,
                            AUTOPAS_MPI_STATUS_IGNORE);

  std::vector<Particle> deserializedParticles{};
  if (not serializedParticles.empty()) {
    Serialization::deserializeParticles(serializedParticles, deserializedParticles);
  }
  return deserializedParticles;
}
