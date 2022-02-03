/**
 * @file Particle.h
 * @author F. Gratl
 * @date 28.06.21
 */

#pragma once

#include <autopas/particles/Particle.h>

#include <array>
#include <ostream>
/**
 * Class describing an arbitrary debris object used for the n-body simulation.
 * Based on the Particle class of AutoPas which used 64bit precision.
 */
class Particle final : public autopas::ParticleFP64 {
 public:
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
   * @param pos
   * @param v
   * @param debrisId
   * @param activityState
   */
  Particle(std::array<double, 3> pos,
           std::array<double, 3> v,
           size_t debrisId,
           ActivityState activityState,
           double mass,
           double radius,
           double coefficientOfDrag)
      : autopas::ParticleFP64(pos, v, debrisId),
        aom(M_PI * radius * radius * 1e-6 / mass),  // convert m^2 -> km^2
        mass(mass),
        radius(radius),
        bc_inv(coefficientOfDrag * aom * 1e6),  // convert km^2 -> m^2
        activityState(activityState) {}

  /**
   * Destructor.
   */
  ~Particle() final = default;

  /**
   * Enums used as ids for accessing and creating a dynamically sized SoA.
   */
  enum AttributeNames : int { ptr, id, posX, posY, posZ, forceX, forceY, forceZ, ownershipState };

  /**
   * The type for the SoA storage.
   */
  using SoAArraysType = typename autopas::utils::SoAType<Particle *,
                                                         size_t /*id*/,
                                                         double /*x*/,
                                                         double /*y*/,
                                                         double /*z*/,
                                                         double /*fx*/,
                                                         double /*fy*/,
                                                         double /*fz*/,
                                                         autopas::OwnershipState /*ownershipState*/>::Type;

  /**
   * Getter, which allows access to an attribute using the corresponding attribute name (defined in AttributeNames).
   * @tparam attribute Attribute name.
   * @return Value of the requested attribute.
   */
  template <AttributeNames attribute>
  [[nodiscard]] constexpr auto get() {
    if constexpr (attribute == AttributeNames::ptr) {
      return this;
    } else if constexpr (attribute == AttributeNames::id) {
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
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      return this->_ownershipState;
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
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      this->_ownershipState = value;
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
   * Calculates speed of the debris
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
  ActivityState getActivityState() const;

  /**
   * Setter for the activityState.
   * @param activityState
   */
  void setActivityState(ActivityState activityState);

  /**
   * Getter for mass.
   * @return
   */
  double getMass() const;

  /**
   * Setter for mass.
   * @param mass
   */
  void setMass(double mass);

  /**
   * Getter for radius.
   * @return
   */
  double getRadius() const;

  /**
   * Setter for radius.
   * @param radius
   */
  void setRadius(double radius);

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

 private:
  /**
   *  3D vector representation of the debris acceleration at the last time step. [km/s^2]
   */
  std::array<double, 3> acc_t0{};
  /**
   * 3D vector representation of the debris acceleration at the current time step. [km/s^2]
   */
  std::array<double, 3> acc_t1{};
  /**
   * Area to mass relation [km^2/kg].
   * @note Used in the Propagator in SRPComponent::apply().
   */
  double aom;
  /**
   * Mass of the object [kg].
   */
  double mass;
  /**
   * Radius of the object when approximating it as a ball [m].
   */
  double radius;
  /**
   * c_D*A/m is the inverse of the ballistic coefficient. [kg m^2]
   * @note Used in the Propagator in Acceleration::DragComponent::apply().
   */
  double bc_inv;

  /**
   * If a particle can actively influence its orbit. See ActivityState.
   */
  ActivityState activityState;
};

/**
 * Input stream operator for the activity state.
 * @note Needed for the CSVReader.
 * @param s
 * @param state
 * @return
 */
std::istream &operator>>(std::istream &s, Particle::ActivityState &state);