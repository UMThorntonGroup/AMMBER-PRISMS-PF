#pragma once

#include <deal.II/base/tensor.h>

template <typename T>
struct is_dealii_tensor
{
  static constexpr bool value = false;
};

template <int rank, int dim, typename T>
struct is_dealii_tensor<dealii::Tensor<rank, dim, T>>
{
  static constexpr bool value = true;
};

template <int dim, typename T>
struct is_dealii_tensor<dealii::Point<dim, T>>
{
  static constexpr bool value = true;
};

template <int rank, int dim, typename T>
struct is_dealii_tensor<dealii::SymmetricTensor<rank, dim, T>>
{
  static constexpr bool value = true;
};

/**
 * @brief Class for holding fields and their spatial derivatives OR variations in their
 * strong form
 * @details Variations are in strong form: [val + div(grad)]
 * @tparam dim The dimension of the problem
 */
template <typename number, unsigned int dim>
struct Dual
{
  using scalarValue = number;
  using scalarGrad  = dealii::Tensor<1, dim, number>;

  scalarValue val  = number(0.);
  scalarGrad  grad = {};

  Dual(const scalarValue &_val, const scalarGrad &_grad)
    : val(_val)
    , grad(_grad)
  {}

  Dual() = default;

  template <typename T>
  explicit Dual(const T &initial_value)
    : val(initial_value)
    , grad()
  {}

  Dual<number, dim>
  operator+(const Dual<number, dim> &other) const
  {
    return {val + other.val, grad + other.grad};
  }

  Dual<number, dim>
  operator+(const scalarValue &constant) const
  {
    return {val + constant, grad};
  }

  Dual<number, dim>
  operator+(const double &constant) const
  {
    return {val + constant, grad};
  }

  Dual<number, dim>
  operator+() const
  {
    return *this;
  }

  template <typename other_number>
  Dual<number, dim> &
  operator+=(const other_number &other)
  {
    *this = *this + other;
    return *this;
  }

  Dual<number, dim>
  operator-(const Dual<number, dim> &other) const
  {
    return {val - other.val, grad - other.grad};
  }

  Dual<number, dim>
  operator-() const
  {
    return {-val, -grad};
  }

  Dual<number, dim>
  operator-(const scalarValue &constant) const
  {
    return {val - constant, grad};
  }

  Dual<number, dim>
  operator-(const double &constant) const
  {
    return {val - constant, grad};
  }

  template <typename other_number>
  Dual<number, dim> &
  operator-=(const other_number &other)
  {
    *this = *this - other;
    return *this;
  }

  Dual<number, dim>
  operator*(const scalarValue &constant) const
  {
    return {val * constant, grad * constant};
  }

  Dual<number, dim>
  operator*(const double &constant) const
  {
    return {val * constant, grad * constant};
  }

  Dual<number, dim>
  operator*(const Dual<number, dim> &other) const
  {
    return Dual<number, dim>(val * other.val, (val * other.grad) + (grad * other.val));
  }

  template <typename other_number, typename = std::enable_if<!is_dealii_tensor<other_number>::value>::type>
  Dual<number, dim> &
  operator*=(const other_number &other)
  {
    *this = *this * other;
    return *this;
  }

  Dual<number, dim>
  operator/(const scalarValue &constant) const
  {
    return Dual<number, dim> {val / constant, grad / constant};
  }

  Dual<number, dim>
  operator/(const double &constant) const
  {
    return Dual<number, dim> {val / constant, grad / constant};
  }

  Dual<number, dim>
  operator/(const Dual<number, dim> &other) const
  {
    return Dual<number, dim> {val / other.val,
                              (grad * other.val - val * other.grad) / (other.val * other.val)};
  }

  template <typename other_number>
  Dual<number, dim> &
  operator/=(const other_number &other)
  {
    *this = *this / other;
    return *this;
  }
};

template <typename number, unsigned int dim, typename other_number>
Dual<number, dim>
operator+(const other_number &other, const Dual<number, dim> &field)
{
  return field + other;
}

template <typename number, unsigned int dim, typename other_number>
Dual<number, dim>
operator-(const other_number &other, const Dual<number, dim> &field)
{
  return -field + other;
}

template <typename number,
          unsigned int dim,
          typename other_number,
          typename = std::enable_if<!is_dealii_tensor<other_number>::value>::type>
Dual<number, dim>
operator*(const other_number &other, const Dual<number, dim> &field)
{
  return field * other;
}

template <typename number, unsigned int dim, typename other_number>
Dual<number, dim>
operator/(const other_number &other, const Dual<number, dim> &field)
{
  Dual<number, dim> one = Dual<number, dim> {number(1.0), {}};
  return (one / field) * other;
}

template <typename number, unsigned int dim>
Dual<number, dim>
sqrt(const Dual<number, dim> &field)
{
  using std::sqrt;
  number sqrt_val = sqrt(field.val);
  return Dual<number, dim> {sqrt_val, field.grad / (2.0 * sqrt_val)};
}

template <typename number, unsigned int dim, typename OStream>
OStream &
operator<<(OStream &stream, const Dual<number, dim> &field)
{
  stream << "Value: " << field.val << ", Gradient: " << field.grad;
  return stream;
}

template <typename number, unsigned int dim>
struct Variation
{
  using scalarValue = number;
  using scalarGrad  = dealii::Tensor<1, dim, number>;

  scalarValue val = number(0.);
  scalarGrad  vec = {};

  Variation(const scalarValue &_val, const scalarGrad &_vec)
    : val(_val)
    , vec(_vec)
  {}

  Variation() = default;

  template <typename T>
  explicit Variation(const T &initial_value)
    : val(initial_value)
    , vec()
  {}

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

template <typename number, unsigned int dim, typename OStream>
OStream &
operator<<(OStream &stream, const Variation<number, dim> &field)
{
  stream << "Value: " << field.val << ", Vector: " << field.vec;
  return stream;
}

namespace dealii
{
  template <typename number, unsigned int dim>
  struct EnableIfScalar<Dual<number, dim>>
  {
    using type = Dual<number, dim>;
  };
} // namespace dealii
