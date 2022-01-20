/**
 * @file Particle.cpp
 * @author F. Gratl
 * @date 28.06.21
 */

#include "Particle.h"

#include "satellitePropagator/utils/MathUtils.h"

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

std::ostream &operator<<(std::ostream &os, const Particle &particle) {
  // clang-format off
  os << particle.toString()
     << "\nMass          : " << particle.getMass()
     << "\nRadius        : " << particle.getRadius()
     << "\nActivityState : " << static_cast<int>(particle.getActivityState());
  // clang-format on
  return os;
}

bool Particle::operator==(const Particle &rhs) const {
  // clang-format off
  return static_cast<const autopas::ParticleBase<double, unsigned long> &>(*this)
      == static_cast<const autopas::ParticleBase<double, unsigned long> &>(rhs) and
      acc_t0 == rhs.acc_t0 and
      acc_t1 == rhs.acc_t1 and
      aom == rhs.aom and
      mass == rhs.mass and
      radius == rhs.radius and
      bc_inv == rhs.bc_inv and
      activityState == rhs.activityState;
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
