#pragma once

#include <deal.II/base/tensor.h>

#include "dual.h"

template <typename number, unsigned int dim>
struct Variation
{
  using scalarValue = number;
  using scalarGrad  = dealii::Tensor<1, dim, number>;

  scalarValue val = number(0.);
  scalarGrad  vec = {};

  template <typename _valuetype = scalarValue, typename _gradtype = scalarGrad>
  Variation(const _valuetype &_val = number(0.), const _gradtype &_vec = {})
    : val(_val)
    , vec(_vec)
  {}

  // Variation(const scalarValue &_val, const scalarGrad &_vec)
  //   : val(_val)
  //   , vec(_vec)
  //{}

  // explicit Variation(const int &initial_value = 0)
  //   : val(number(double(initial_value)))
  //   , vec()
  //{}

  Variation<number, dim>
  operator+(const Variation<number, dim> &other) const
  {
    return Variation<number, dim> {val + other.val, vec + other.vec};
  }

  Variation<number, dim>
  operator+(const scalarValue &scalar) const
  {
    return Variation<number, dim> {val + scalar, vec};
  }

  Variation<number, dim>
  operator+(const double &scalar) const
  {
    return Variation<number, dim> {val + scalar, vec};
  }

  Variation<number, dim>
  operator+(const scalarGrad &vector) const
  {
    return Variation<number, dim> {val, vec + vector};
  }

  Variation<number, dim>
  operator+() const
  {
    return *this;
  }

  template <typename other_number>
  Variation<number, dim> &
  operator+=(const other_number &other)
  {
    *this = *this + other;
    return *this;
  }

  Variation<number, dim>
  operator-(const Variation<number, dim> &other) const
  {
    return Variation<number, dim> {val - other.val, vec - other.vec};
  }

  Variation<number, dim>
  operator-(const scalarValue &scalar) const
  {
    return Variation<number, dim> {val - scalar, vec};
  }

  Variation<number, dim>
  operator-(const double &scalar) const
  {
    return Variation<number, dim> {val - scalar, vec};
  }

  Variation<number, dim>
  operator-(const scalarGrad &vector) const
  {
    return Variation<number, dim> {val, vec - vector};
  }

  Variation<number, dim>
  operator-() const
  {
    return Variation<number, dim> {-val, -vec};
  }

  template <typename other_number>
  Variation<number, dim> &
  operator-=(const other_number &other)
  {
    *this = *this - other;
    return *this;
  }

  Variation<number, dim>
  operator*(const scalarValue &constant_scalar) const
  {
    return Variation<number, dim> {val * constant_scalar, vec * constant_scalar};
  }

  Variation<number, dim>
  operator*(const double &constant_scalar) const
  {
    return Variation<number, dim> {val * constant_scalar, vec * constant_scalar};
  }

  template <typename other_number>
  Variation<number, dim> &
  operator*=(const other_number &other)
  {
    *this = *this * other;
    return *this;
  }

  Variation<number, dim>
  operator/(const scalarValue &constant_scalar) const
  {
    return Variation<number, dim> {val / constant_scalar, vec / constant_scalar};
  }

  Variation<number, dim>
  operator/(const double &constant_scalar) const
  {
    return Variation<number, dim> {val / constant_scalar, vec / constant_scalar};
  }

  template <typename other_number>
  Variation<number, dim> &
  operator/=(const other_number &other)
  {
    *this = *this / other;
    return *this;
  }
};

template <typename number, unsigned int dim, typename other_number>
Variation<number, dim>
operator+(const other_number &other, const Variation<number, dim> &variation)
{
  return variation + other;
}

template <typename number, unsigned int dim, typename other_number>
Variation<number, dim>
operator-(const other_number &other, const Variation<number, dim> &variation)
{
  return -variation + other;
}

template <typename number, unsigned int dim, typename other_number>
Variation<number, dim>
operator*(const other_number &other, const Variation<number, dim> &variation)
{
  return variation * other;
}

/**
 * @brief Multiply a field and a variation
 * @param field The field to multiply
 * @param variation The variation to multiply
 */
template <typename number, unsigned int dim>
Variation<number, dim>
operator*(const Dual<number, dim> &field, const Variation<number, dim> &variation)
{
  return Variation<number, dim> {field.val * variation.val - field.grad * variation.vec,
                                 field.val * variation.vec};
}

/**
 * @brief Multiply a field and a variation
 * @param field The field to multiply
 * @param variation The variation to multiply
 */
template <typename number, unsigned int dim>
Variation<number, dim>
operator*(const Variation<number, dim> &variation, const Dual<number, dim> &field)
{
  return field * variation;
}