/**
 * @file Particle.h
 * @author F. Gratl
 * @date 28.06.21
 */

#pragma once

#include <autopas/particles/Particle.h>

#include <array>
#include <ostream>
#include <utility>

namespace LADDS {

/**
 * Class describing an arbitrary debris object used for the n-body simulation.
 * Based on the Particle class of AutoPas which used 64bit precision.
 */
class Particle final : public autopas::ParticleFP64 {
 public:
  /**
   * calculate the inverse ballistic coefficient from bstar if available,
   * otherwise fall back to radius and mass.
   *
   * @param bstar [1/R_EARTH]
   * @param radius [m]
   * @param mass [kg]
   * @param coefficientOfDrag
   */
  static double calculateBcInv(double bstar, double mass, double radius, double coefficientOfDrag);

  /**
   * Describes how active a particle behaves. This is relevant for the propagator to determine which forces are applied.
   */
  enum class ActivityState : int {
    /**
     * Simply float around space and is subject to all external influences.
     */
    passive,

    /**
     * Can actively change its trajectory to avoid collisions.
     */
    evasive,

    /**
     * Same as evasive but additionally can actively maneuver to preserve their orbit.
     * This is modelled only only applying Keplerian forces.
     */
    evasivePreserving,
  };

  /**
   * Constructor
   * @param pos [km]
   * @param v [km/s]
   * @param debrisId
   * @param activityState
   */
  Particle(std::array<double, 3> pos,
           std::array<double, 3> v,
           size_t debrisId,
           std::string identifier,
           ActivityState activityState,
           double mass,
           double radius,
           double bcInv,
           size_t parentIdentifier)
      : autopas::ParticleFP64(pos, v, debrisId),
        _aom(M_PI * radius * radius * 1e-6 / mass),  // convert m^2 -> km^2
        _mass(mass),
        _radius(radius),
        _bc_inv(bcInv),
        _activityState(activityState),
        _identifier(std::move(identifier)),
        _parentIdentifier(parentIdentifier) {
    // to make serialization easier make sure COSPAR IDs are always padded with spaces to full length of 11
    _identifier.insert(_identifier.end(), 11 - _identifier.size(), ' ');
    // FIXME: remove this check as soon as we are confident that it works
    if (_identifier.length() > 11) {
      std::cerr << "WARNING identifier.length() == " << _identifier.length() << " |" << _identifier << "|" << std::endl;
    }
  }

  /**
   * Destructor.
   */
  ~Particle() final = default;

 private:
  /**
   *  3D vector representation of the debris acceleration at the last time step. [km/s^2]
   */
  std::array<double, 3> _acc_t0{};
  /**
   * 3D vector representation of the debris acceleration at the current time step. [km/s^2]
   */
  std::array<double, 3> _acc_t1{};
  /**
   * Area to _mass relation [km^2/kg].
   * @note Used in the Propagator in SRPComponent::apply().
   */
  double _aom;
  /**
   * Mass of the object [kg].
   */
  double _mass;
  /**
   * Radius of the object when approximating it as a ball [m].
   */
  double _radius;
  /**
   * c_D*A/m is the inverse of the ballistic coefficient. [kg m^2]
   * @note Used in the Propagator in Acceleration::DragComponent::apply().
   */
  double _bc_inv;

  /**
   * If a particle can actively influence its orbit. See ActivityState.
   */
  ActivityState _activityState;

  /**
   * Unique string _identifier to relate objects to catalogue objects
   */
  std::string _identifier;

  /**
   * Unique identifier of the particle whose collision created this one
   */
  size_t _parentIdentifier;

 public:
  /**
   * Enums used as ids for accessing and creating a dynamically sized SoA.
   */
  enum AttributeNames : int {
    ptr,
    id,
    posX,
    posY,
    posZ,
    forceX,
    forceY,
    forceZ,
    velocityX,
    velocityY,
    velocityZ,
    ownershipState,
    acc_t0X,
    acc_t0Y,
    acc_t0Z,
    acc_t1X,
    acc_t1Y,
    acc_t1Z,
    aom,
    mass,
    radius,
    bc_inv,
    activityState,
    identifier,
    parentIdentifier
  };

  /**
   * The type for the SoA storage.
   */
  using SoAArraysType = typename autopas::utils::SoAType<Particle *,
                                                         decltype(_id) /*id*/,
                                                         std::remove_reference_t<decltype(_r[0])> /*x*/,
                                                         std::remove_reference_t<decltype(_r[1])> /*y*/,
                                                         std::remove_reference_t<decltype(_r[2])> /*z*/,
                                                         std::remove_reference_t<decltype(_f[0])> /*fx*/,
                                                         std::remove_reference_t<decltype(_f[1])> /*fy*/,
                                                         std::remove_reference_t<decltype(_f[2])> /*fz*/,
                                                         std::remove_reference_t<decltype(_v[0])> /*vx*/,
                                                         std::remove_reference_t<decltype(_v[1])> /*vy*/,
                                                         std::remove_reference_t<decltype(_v[2])> /*vz*/,
                                                         decltype(_ownershipState) /*ownershipState*/,
                                                         std::remove_reference_t<decltype(_acc_t0[0])> /*acc_t0X*/,
                                                         std::remove_reference_t<decltype(_acc_t0[1])> /*acc_t0Y*/,
                                                         std::remove_reference_t<decltype(_acc_t0[2])> /*acc_t0Z*/,
                                                         std::remove_reference_t<decltype(_acc_t1[0])> /*acc_t1X*/,
                                                         std::remove_reference_t<decltype(_acc_t1[1])> /*acc_t1Y*/,
                                                         std::remove_reference_t<decltype(_acc_t1[2])> /*acc_t1Z*/,
                                                         decltype(_aom) /*_aom*/,
                                                         decltype(_mass) /*_mass*/,
                                                         decltype(_radius) /*_radius*/,
                                                         decltype(_bc_inv) /*_bc_inv*/,
                                                         decltype(_activityState) /*_activityState*/,
                                                         decltype(_identifier) /*_identifier*/,
                                                         decltype(_parentIdentifier) /*_parentIdentifier*/
                                                         >::Type;

  /**
   * Non-const getter for the pointer of this object.
   * @tparam attribute Attribute name.
   * @return this.
   */
  template <AttributeNames attribute, std::enable_if_t<attribute == AttributeNames::ptr, bool> = true>
  [[nodiscard]] constexpr auto get() {
    return this;
  }

  /**
   * Getter, which allows access to an attribute using the corresponding attribute name (defined in AttributeNames).
   * @tparam attribute Attribute name.
   * @return Value of the requested attribute.
   */
  template <AttributeNames attribute, std::enable_if_t<attribute != AttributeNames::ptr, bool> = true>
  [[nodiscard]] constexpr auto get() const {
    if constexpr (attribute == AttributeNames::id) {
      return getID();
    } else if constexpr (attribute == AttributeNames::posX) {
      return getR()[0];
    } else if constexpr (attribute == AttributeNames::posY) {
      return getR()[1];
    } else if constexpr (attribute == AttributeNames::posZ) {
      return getR()[2];
    } else if constexpr (attribute == AttributeNames::forceX) {
      return getF()[0];
    } else if constexpr (attribute == AttributeNames::forceY) {
      return getF()[1];
    } else if constexpr (attribute == AttributeNames::forceZ) {
      return getF()[2];
    } else if constexpr (attribute == AttributeNames::velocityX) {
      return getVelocity()[0];
    } else if constexpr (attribute == AttributeNames::velocityY) {
      return getVelocity()[1];
    } else if constexpr (attribute == AttributeNames::velocityZ) {
      return getVelocity()[2];
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      return getOwnershipState();
    } else if constexpr (attribute == AttributeNames::acc_t0X) {
      return getAccT0()[0];
    } else if constexpr (attribute == AttributeNames::acc_t0Y) {
      return getAccT0()[1];
    } else if constexpr (attribute == AttributeNames::acc_t0Z) {
      return getAccT0()[2];
    } else if constexpr (attribute == AttributeNames::acc_t1X) {
      return getAccT1()[0];
    } else if constexpr (attribute == AttributeNames::acc_t1Y) {
      return getAccT1()[1];
    } else if constexpr (attribute == AttributeNames::acc_t1Z) {
      return getAccT1()[2];
    } else if constexpr (attribute == AttributeNames::aom) {
      return getAom();
    } else if constexpr (attribute == AttributeNames::mass) {
      return getMass();
    } else if constexpr (attribute == AttributeNames::radius) {
      return getRadius();
    } else if constexpr (attribute == AttributeNames::bc_inv) {
      return getBcInv();
    } else if constexpr (attribute == AttributeNames::activityState) {
      return getActivityState();
    } else if constexpr (attribute == AttributeNames::identifier) {
      return getIdentifier();
    } else {
      autopas::utils::ExceptionHandler::exception("Particle::get() unknown attribute {}", attribute);
    }
  }

  /**
   * Setter, which allows set an attribute using the corresponding attribute name (defined in AttributeNames).
   * @tparam attribute Attribute name.
   * @param value New value of the requested attribute.
   */
  template <AttributeNames attribute>
  constexpr void set(typename std::tuple_element<attribute, SoAArraysType>::type::value_type value) {
    if constexpr (attribute == AttributeNames::id) {
      setID(value);
    } else if constexpr (attribute == AttributeNames::posX) {
      _r[0] = value;
    } else if constexpr (attribute == AttributeNames::posY) {
      _r[1] = value;
    } else if constexpr (attribute == AttributeNames::posZ) {
      _r[2] = value;
    } else if constexpr (attribute == AttributeNames::forceX) {
      _f[0] = value;
    } else if constexpr (attribute == AttributeNames::forceY) {
      _f[1] = value;
    } else if constexpr (attribute == AttributeNames::forceZ) {
      _f[2] = value;
    } else if constexpr (attribute == AttributeNames::velocityX) {
      _v[0] = value;
    } else if constexpr (attribute == AttributeNames::velocityY) {
      _v[1] = value;
    } else if constexpr (attribute == AttributeNames::velocityZ) {
      _v[2] = value;
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      this->_ownershipState = value;
    } else if constexpr (attribute == AttributeNames::acc_t0X) {
      _acc_t0[0] = value;
    } else if constexpr (attribute == AttributeNames::acc_t0Y) {
      _acc_t0[1] = value;
    } else if constexpr (attribute == AttributeNames::acc_t0Z) {
      _acc_t0[2] = value;
    } else if constexpr (attribute == AttributeNames::acc_t1X) {
      _acc_t1[0] = value;
    } else if constexpr (attribute == AttributeNames::acc_t1Y) {
      _acc_t1[1] = value;
    } else if constexpr (attribute == AttributeNames::acc_t1Z) {
      _acc_t1[2] = value;
    } else if constexpr (attribute == AttributeNames::aom) {
      _aom = value;
    } else if constexpr (attribute == AttributeNames::mass) {
      _mass = value;
    } else if constexpr (attribute == AttributeNames::radius) {
      _radius = value;
    } else if constexpr (attribute == AttributeNames::bc_inv) {
      _bc_inv = value;
    } else if constexpr (attribute == AttributeNames::activityState) {
      _activityState = value;
    } else if constexpr (attribute == AttributeNames::identifier) {
      _identifier = value;
    } else if constexpr (attribute == AttributeNames::parentIdentifier) {
      _parentIdentifier = value;
    } else {
      autopas::utils::ExceptionHandler::exception("Particle::set() unknown attribute {}", attribute);
    }
  }

  /**
   * Calculates distance from the origin of the coordinate frame [km]
   *
   * @return Euclidean norm of the #position vector
   */
  [[nodiscard]] double getHeight() const;

  /**
   * Calculates speed of the debris [km/s]
   *
   * @return Euclidean norm of the #velocity vector
   */
  [[nodiscard]] double getSpeed() const;

  /**
   * Calculates the euclidean norm of the #acc_t0
   *
   * @return Calculates the euclidean norm of the #acc_t0
   */
  [[nodiscard]] double getAccT0Norm() const;

  /**
   * Calculates the euclidean norm of the #acc_t1
   *
   * @return Calculates the euclidean norm of the #acc_t1
   */
  [[nodiscard]] double getAccT1Norm() const;

  /**
   * Getter function for #position vector [km]
   *
   * @return 3D vector representation of the debris #position
   */
  [[nodiscard]] const std::array<double, 3> &getPosition() const;

  /**
   * Setter function for #position vector [km]
   *
   * @param position 3D vector representation of the debris #position
   */
  void setPosition(const std::array<double, 3> &position);

  /**
   * Getter function for #velocity vector [km/s]
   *
   * @return 3D vector representation of the debris #velocity
   */
  [[nodiscard]] const std::array<double, 3> &getVelocity() const;

  /**
   * Setter function for #velocity vector
   *
   * @param velocity 3D vector representation of the debris #velocity
   */
  void setVelocity(const std::array<double, 3> &velocity);

  /**
   * Getter function for #acc_t0 vector
   *
   * @return 3D vector representation of the debris #acc_t0
   */
  [[nodiscard]] const std::array<double, 3> &getAccT0() const;

  /**
   * Setter function for #acc_t0 vector
   *
   * @param accT0 3D vector representation of the debris #acc_t0
   */
  void setAccT0(const std::array<double, 3> &accT0);

  /**
   * Getter function for #acc_t1 vector
   *
   * @return 3D vector representation of the debris #acc_t1
   */
  [[nodiscard]] const std::array<double, 3> &getAccT1() const;

  /**
   * Setter function for #acc_t1 vector
   *
   * @param accT1 3D vector representation of the debris #acc_t1
   */
  void setAccT1(const std::array<double, 3> &accT1);

  /**
   * Getter function for #aom
   *
   * @return value of #aom
   */
  [[nodiscard]] double getAom() const;

  /**
   * Setter function for #aom
   *
   * @param aom New value #aom
   */
  void setAom(double aom);

  /**
   * Getter function for #bc_inv
   *
   * @return Value of #bc_inv
   */
  [[nodiscard]] double getBcInv() const;

  /**
   * Setter function for #bc_inv
   *
   * @param bcInv New value of#bc_inv
   */
  void setBcInv(double bcInv);

  /**
   * Getter for the activityState.
   * @return
   */
  [[nodiscard]] ActivityState getActivityState() const;

  /**
   * Setter for the activityState.
   * @param activityState
   */
  void setActivityState(ActivityState activityState);

  /**
   * Getter for mass.
   * @return
   */
  [[nodiscard]] double getMass() const;

  /**
   * Setter for mass.
   * @param mass
   */
  void setMass(double mass);

  /**
   * Getter for radius.
   * @return
   */
  [[nodiscard]] double getRadius() const;

  /**
   * Setter for radius.
   * @param radius
   */
  void setRadius(double radius);

  /**
   * Getter for identifier.
   * @return
   */
  [[nodiscard]] const std::string &getIdentifier() const;

  /**
   * Setter for identifier.
   * @param
   */
  void setIdentifier(const std::string &identifier);

  /**
   * Getter for identifier.
   * @return
   */
  [[nodiscard]] const size_t &getParentIdentifier() const;

  /**
   * Setter for identifier.
   * @param
   */
  void setParentIdentifier(const size_t &parentIdentifier);

  /**
   * Stream operator
   * @param os
   * @param particle
   * @return
   */
  friend std::ostream &operator<<(std::ostream &os, const Particle &particle);

  /**
   * Equality operator
   * @param rhs
   * @return
   */
  bool operator==(const Particle &rhs) const;

  /**
   * Inequality operator
   * @param rhs
   * @return
   */
  bool operator!=(const Particle &rhs) const;
};

/**
 * Input stream operator for the activity state.
 * @note Needed for the CSVReader.
 * @param s
 * @param state
 * @return
 */
std::istream &operator>>(std::istream &s, Particle::ActivityState &state);

}  // namespace LADDS