/**
 * @file Serialization.h
 * @author F. Gratl
 * @date 31.05.22
 */
#pragma once

#include <string>
#include <utility>
#include <vector>

#include "ladds/particle/Particle.h"
/**
 * Provides tools to de-/serialize particles.
 */
namespace LADDS::ParticleSerializationTools {
/**
 * Serializes a particle and appends it to the serializedParticles container.
 * @param particle The particle which will be serialized.
 * @param serializedParticles The container to wich the serialized particle will be appended.
 */
void serializeParticle(const Particle &particle, std::vector<char> &serializedParticles);

/**
 * Deserializes a serialized particle.
 * @param particleData A pointer to the serialized particle data.
 * @param particle The particle to which the desierialized attributes will be applied.
 */
void deserializeParticle(char *particleData, Particle &particle);

/**
 * Deserializes a container of serialized particles.
 * @param particlesData A pointer to the serialized particle data.
 * @param particles The particle container to which to append the deserialized particles to.
 */
void deserializeParticles(std::vector<char> &particlesData, std::vector<Particle> &particles);
}  // namespace LADDS::ParticleSerializationTools
