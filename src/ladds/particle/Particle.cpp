/**
 * @file Particle.cpp
 * @author F. Gratl
 * @date 28.06.21
 */
#include "Particle.h"

#include <satellitePropagator/physics/Constants.h>
#include <satellitePropagator/utils/MathUtils.h>

namespace LADDS {

double Particle::calculateBcInv(double bstar, double mass, double radius, double coefficientOfDrag) {
  if (std::isnan(bstar) or bstar == 0.) {
    // either via c_D
    const auto area = M_PI * radius * radius;  // m^2
    return coefficientOfDrag * area / mass;    // m^2/kg
  } else {
    // or via bstar
    // @note see https://en.wikipedia.org/wiki/BSTAR
    // B* == p_0 * c_D * A / (2 m) == bc_inv * p_0 / 2
    // Thus we factor out p0 (=0.1570), convert to m (original bstar is in Earth radii) and the factor 2 comes
    // from the propagator still having a 1/2 multiplication term included which is already accounted for in BSTAR
    constexpr auto factor = 2.0 * Physics::R_EARTH / 0.1570;
    return factor * bstar;
  }
}

double Particle::getHeight() const {
  return MathUtils::euclideanNorm(getPosition());
}

double Particle::getSpeed() const {
  return MathUtils::euclideanNorm(getVelocity());
}

double Particle::getAccT0Norm() const {
  return MathUtils::euclideanNorm(acc_t0);
}

double Particle::getAccT1Norm() const {
  return MathUtils::euclideanNorm(acc_t1);
}

const std::array<double, 3> &Particle::getPosition() const {
  return _r;
}

void Particle::setPosition(const std::array<double, 3> &position) {
  Particle::_r = position;
}

const std::array<double, 3> &Particle::getVelocity() const {
  return _v;
}

void Particle::setVelocity(const std::array<double, 3> &velocity) {
  Particle::_v = velocity;
}

const std::array<double, 3> &Particle::getAccT0() const {
  return acc_t0;
}

void Particle::setAccT0(const std::array<double, 3> &accT0) {
  acc_t0 = accT0;
}

const std::array<double, 3> &Particle::getAccT1() const {
  return acc_t1;
}

void Particle::setAccT1(const std::array<double, 3> &accT1) {
  acc_t1 = accT1;
}

double Particle::getAom() const {
  return aom;
}

void Particle::setAom(const double aom) {
  Particle::aom = aom;
}

double Particle::getBcInv() const {
  return bc_inv;
}

void Particle::setBcInv(const double bcInv) {
  bc_inv = bcInv;
}

Particle::ActivityState Particle::getActivityState() const {
  return activityState;
}

void Particle::setActivityState(Particle::ActivityState activityState) {
  Particle::activityState = activityState;
}

const std::string &Particle::getIdentifier() const {
  return identifier;
}

void Particle::setIdentifier(const std::string &identifier) {
  this->identifier = identifier;
}

std::ostream &operator<<(std::ostream &os, const Particle &particle) {
  // clang-format off
  os << particle.toString()
     << "\nIdentifier    : " << particle.getIdentifier()
     << "\nMass          : " << particle.getMass()
     << "\nRadius        : " << particle.getRadius()
     << "\nBcInv         : " << particle.getBcInv()
     << "\nActivityState : " << static_cast<int>(particle.getActivityState());
  // clang-format on
  return os;
}

std::istream &operator>>(std::istream &s, Particle::ActivityState &state) {
  std::string str;
  s >> str;
  if (str == "passive") {
    state = Particle::ActivityState::passive;
  } else if (str == "evasive") {
    state = Particle::ActivityState::evasive;
  } else if (str == "evasivePreserving") {
    state = Particle::ActivityState::evasivePreserving;
  }
  return s;
}

bool Particle::operator==(const Particle &rhs) const {
  // clang-format off
  return static_cast<const autopas::ParticleBase<double, unsigned long> &>(*this)
      == static_cast<const autopas::ParticleBase<double, unsigned long> &>(rhs) and
      acc_t0 == rhs.acc_t0 and
      acc_t1 == rhs.acc_t1 and
      // the following values can only be compared upon float32 bit precision since we write them as such to HDF5
      (std::abs(aom - rhs.aom) < 1e-7) and
      (std::abs(mass - rhs.mass) < 1e-7) and
      (std::abs(radius - rhs.radius) < 1e-7) and
      (std::abs(bc_inv - rhs.bc_inv) < 1e-7) and
      activityState == rhs.activityState and
      identifier == rhs.identifier;
  // clang-format on
}
bool Particle::operator!=(const Particle &rhs) const {
  return !(rhs == *this);
}

double Particle::getMass() const {
  return mass;
}

void Particle::setMass(double mass) {
  Particle::mass = mass;
}

double Particle::getRadius() const {
  return radius;
}

void Particle::setRadius(double radius) {
  Particle::radius = radius;
}

}  // namespace LADDS
