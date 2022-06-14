/**
 * @file Serialization.cpp
 * @author F. Gratl
 * @date 31.05.22
 */

#include "Serialization.h"

#include <tuple>

namespace LADDS {
namespace {

/**
 * Stores the AttributeNames of the attributes of ParticleType which have to be communicated using MPI.
 */
constexpr std::array<typename Particle::AttributeNames, 14> Attributes = {Particle::AttributeNames::id,
                                                                          Particle::AttributeNames::posX,
                                                                          Particle::AttributeNames::posY,
                                                                          Particle::AttributeNames::posZ,
                                                                          // Particle::AttributeNames::forceX,
                                                                          // Particle::AttributeNames::forceY,
                                                                          // Particle::AttributeNames::forceZ,
                                                                          Particle::AttributeNames::velocityX,
                                                                          Particle::AttributeNames::velocityY,
                                                                          Particle::AttributeNames::velocityZ,
                                                                          Particle::AttributeNames::ownershipState,
                                                                          // Particle::AttributeNames::acc_t0X,
                                                                          // Particle::AttributeNames::acc_t0Y,
                                                                          // Particle::AttributeNames::acc_t0Z,
                                                                          // Particle::AttributeNames::acc_t1X,
                                                                          // Particle::AttributeNames::acc_t1Y,
                                                                          // Particle::AttributeNames::acc_t1Z,
                                                                          Particle::AttributeNames::aom,
                                                                          Particle::AttributeNames::mass,
                                                                          Particle::AttributeNames::radius,
                                                                          Particle::AttributeNames::bc_inv,
                                                                          Particle::AttributeNames::activityState,
                                                                          Particle::AttributeNames::identifier};

/**
 * The combined size in byte of the attributes which need to be communicated using MPI.
 * Not really sure where the 4 extra bytes in the end come from but they are needed.
 */
constexpr size_t AttributesSize =
    1 * sizeof(size_t) /*id*/ + 10 * sizeof(double) + 2 * sizeof(int) /*enums*/ + 12 * sizeof(char) /*identifier*/ + 4;

/**
 * Serializes the attribute defined by I.
 * @param particle: The particle who's attribute needs to be serialized.
 * @param attributeVector: The container in which the serialized attribute will be stored.
 * @param startIndex: The startindex in the container where to store the serialized attribute.
 */
template <size_t I>
void serializeAttribute(const Particle &particle, std::vector<char> &attributeVector, size_t &startIndex) {
  const auto attribute = particle.get<Attributes[I]>();
  constexpr auto index = Attributes[I];
  if constexpr (index == LADDS::Particle::AttributeNames::identifier) {
    const auto sizeOfValue = 12;
    std::memcpy(&attributeVector[startIndex], attribute.data(), sizeOfValue);
    startIndex += sizeOfValue;
  } else {
    const auto sizeOfValue = sizeof(attribute);
    std::memcpy(&attributeVector[startIndex], &attribute, sizeOfValue);
    startIndex += sizeOfValue;
  }
}

/**
 * Deserializes the attribute defined by I.
 * @param attributeVector: The vector containing the data which needs to be deserialized.
 * @param particle: The particle to which the serialized data will be applied.
 * @param startIndex: The start index in the attributeVector of the attribute which needs to be deserialized.
 */
template <size_t I>
void deserializeAttribute(char *&attributeVector, Particle &particle, size_t &startIndex) {
  auto i = I;
  constexpr auto index = Attributes[I];
  auto attribute = particle.get<index>();
  if constexpr (index == LADDS::Particle::AttributeNames::identifier) {
    const auto sizeOfValue = 12;
    std::memcpy(attribute.data(), &attributeVector[startIndex], sizeOfValue);
    startIndex += sizeOfValue;
  } else {
    const auto sizeOfValue = sizeof(attribute);
    std::memcpy(&attribute, &attributeVector[startIndex], sizeOfValue);
    startIndex += sizeOfValue;
  }
  particle.set<Attributes[I]>(attribute);
}

/**
 * The implementation of serializeParticle using the expansion operator.
 * @param particle: The particle which will be serialized.
 * @param serializedParticle: The char array of the particles serialized attributes.
 */
template <size_t... I>
void serializeParticleImpl(const Particle &particle, std::vector<char> &serializedParticle, std::index_sequence<I...>) {
  // Serialize particle attributes
  size_t startIndex = 0;
  std::vector<char> attributesVector(AttributesSize);
  (serializeAttribute<I>(particle, attributesVector, startIndex), ...);

  // Add serialized attributes to serialized particle
  serializedParticle.insert(serializedParticle.end(), attributesVector.begin(), attributesVector.end());
}

/**
 * The implementation fo deserializeParticle using the expansion operator.
 * @param particleData: The particle data which will be deserialized.
 * @param particle: The particle to which the deserialized attributes will be applied.
 */
template <size_t... I>
void deserializeParticleImpl(char *particleData, Particle &particle, std::index_sequence<I...>) {
  size_t startIndex = 0;
  (deserializeAttribute<I>(particleData, particle, startIndex), ...);
}
}  // namespace

namespace Serialization {

void serializeParticles(const std::vector<Particle>::iterator &particlesBegin,
                        const std::vector<Particle>::iterator &particlesEnd,
                        std::vector<char> &serializedParticles) {
  for (auto particleIter = particlesBegin; particleIter < particlesEnd; ++particleIter) {
    serializeParticleImpl(*particleIter, serializedParticles, std::make_index_sequence<Attributes.size()>{});
  }
}

void deserializeParticles(std::vector<char> &particlesData, std::vector<Particle> &particles) {
  particles.reserve(particles.size() + (particlesData.size() / AttributesSize));
  // initialize a dummy particle which will be filled with parsed values
  Particle particle{{0., 0., 0.}, {0., 0., 0.}, 0, "0123-456abc", Particle::ActivityState::passive, 0., 0., 0.};
  for (size_t i = 0; i < particlesData.size(); i += AttributesSize) {
    deserializeParticleImpl(&particlesData[i], particle, std::make_index_sequence<Attributes.size()>{});
    particles.push_back(particle);
  }
}

}  // namespace Serialization
}  // namespace LADDS