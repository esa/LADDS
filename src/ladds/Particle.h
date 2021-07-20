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
  explicit Particle(std::array<double, 3> pos, std::array<double, 3> v, size_t debrisId)
      : autopas::ParticleFP64(pos, v, debrisId) {}

  ~Particle() final = default;

  /**
   * Enums used as ids for accessing and creating a dynamically sized SoA.
   */
  enum AttributeNames : int { ptr, id, posX, posY, posZ, forceX, forceY, forceZ, ownershipState };

  /**
   * The type for the SoA storage.
   */
  using SoAArraysType =
      typename autopas::utils::SoAType<Particle *, size_t /*id*/, double /*x*/, double /*y*/, double /*z*/, double /*fx*/,
                                       double /*fy*/, double /*fz*/, autopas::OwnershipState /*ownershipState*/>::Type;

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
   * @brief Calculates distance from the origin of the coordinate frame
   *
   * @return Euclidean norm of the #position vector
   */
  [[nodiscard]] double getHeight() const;

  /**
   * @brief Calculates speed of the debris
   *
   * @return Euclidean norm of the #velocity vector
   */
  [[nodiscard]] double getSpeed() const;

  /**
   * @brief Calculates the euclidean norm of the #acc_t0
   *
   * @return Calculates the euclidean norm of the #acc_t0
   */
  [[nodiscard]] double getAccT0Norm() const;

  /**
   * @brief Calculates the euclidean norm of the #acc_t1
   *
   * @return Calculates the euclidean norm of the #acc_t1
   */
  [[nodiscard]] double getAccT1Norm() const;

  /**
   * @brief Getter function for #position vector
   *
   * @return 3D vector representation of the debris #position
   */
  [[nodiscard]] const std::array<double, 3>& getPosition() const;
  //std::array<double, 3>& getPosition();

  /**
   * @brief Setter function for #position vector
   *
   * @param position 3D vector representation of the debris #position
   */
  void setPosition(const std::array<double, 3>& position);

  /**
   * @brief Getter function for #velocity vector
   *
   * @return 3D vector representation of the debris #velocity
   */
  [[nodiscard]] const std::array<double, 3>& getVelocity() const;
  //std::array<double, 3>& getVelocity();

  /**
   * @brief Setter function for #velocity vector
   *
   * @param velocity 3D vector representation of the debris #velocity
   */
  void setVelocity(const std::array<double, 3>& velocity);

  /**
   * @brief Getter function for #acc_t0 vector
   *
   * @return 3D vector representation of the debris #acc_t0
   */
  [[nodiscard]] const std::array<double, 3>& getAccT0() const;
  //std::array<double, 3>& getAccT0();

  /**
   * @brief Setter function for #acc_t0 vector
   *
   * @param accT0 3D vector representation of the debris #acc_t0
   */
  void setAccT0(const std::array<double, 3>& accT0);

  /**
   * @brief Getter function for #acc_t1 vector
   *
   * @return 3D vector representation of the debris #acc_t1
   */
  [[nodiscard]] const std::array<double, 3>& getAccT1() const;
  //std::array<double, 3>& getAccT1();

  /**
   * @brief Setter function for #acc_t1 vector
   *
   * @param accT1 3D vector representation of the debris #acc_t1
   */
  void setAccT1(const std::array<double, 3>& accT1);

  /**
   * @brief Getter function for #aom
   *
   * @return value of #aom
   */
  [[nodiscard]] double getAom() const;

  /**
   * @brief Setter function for #aom
   *
   * @param aom New value #aom
   */
  void setAom(double aom);

  /**
   * @brief Getter function for #bc_inv
   *
   * @return Value of #bc_inv
   */
  [[nodiscard]] double getBcInv() const;

  /**
   * @brief Setter function for #bc_inv
   *
   * @param bcInv New value of#bc_inv
   */
  void setBcInv(double bcInv);
 private:

  /**
   *  3D vector representation of the debris acceleration at the last time step.
   */
  std::array<double, 3> acc_t0 {}; 
  /**
   * 3D vector representation of the debris acceleration at the current time step
   */
  std::array<double, 3> acc_t1 {}; 
  /**
   * Area to mass ration.
   */
  double aom = 0; 
  /**
   * C_cA)/m is the inverse of the ballistic coefficient. 
   * Used for Acceleration::DragComponent::apply().
   */
  double bc_inv = 0; 

};
