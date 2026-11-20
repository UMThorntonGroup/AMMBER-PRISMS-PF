#pragma once

#include <ammber/utilities/auto_diff.h>

#include <deal.II/base/tensor.h>

// Step 1: Cholesky decomposition
template <int n, typename number>
dealii::Tensor<2, n, number>
cholesky_decompose(const dealii::SymmetricTensor<2, n, number> &A)
{
  using std::sqrt;

  dealii::Tensor<2, n, number> L;

  for (std::size_t i = 0; i < n; ++i)
    {
      for (std::size_t j = 0; j <= i; ++j)
        {
          number sum = A[i][j];
          for (std::size_t k = 0; k < j; ++k)
            {
              sum -= L[i][k] * L[j][k];
            }

          if (i == j)
            {
              // if (sum <= 0.0)
              //   return false; // Not positive definite
              L[i][j] = sqrt(sum);
            }
          else
            {
              L[i][j] = sum / L[j][j];
            }
        }
    }
  return L;
}

// Step 2: Forward substitution
template <int n, typename number>
dealii::Tensor<1, n, number>
forward_substitution(const dealii::Tensor<2, n, number> &L, const dealii::Tensor<1, n, number> &b)
{
  dealii::Tensor<1, n, number> y;
  for (std::size_t i = 0; i < n; ++i)
    {
      number sum = b[i];
      for (std::size_t j = 0; j < i; ++j)
        {
          sum -= L[i][j] * y[j];
        }
      y[i] = sum / L[i][i];
    }
  return y;
}

// Step 3: Backward substitution (for L^T)
template <int n, typename number>
dealii::Tensor<1, n, number>
backward_substitution(const dealii::Tensor<2, n, number> &L, const dealii::Tensor<1, n, number> &y)
{
  dealii::Tensor<1, n, number> x;
  for (int i = n - 1; i >= 0; --i)
    {
      number sum = y[i];
      for (int j = i + 1; j < n; ++j)
        {
          sum -= L[j][i] * x[j]; // Note: L^T used here
        }
      x[i] = sum / L[i][i];
    }
  return x;
}

// Step 4: Compute A inverse
template <int n, typename number>
dealii::Tensor<2, n, number>
inverse_from_cholesky(const dealii::SymmetricTensor<2, n, number> &A)
{
  if constexpr (n <= 3)
    {
      return dealii::invert(A);
    }
  dealii::Tensor<2, n, number> L = cholesky_decompose(A);

  dealii::Tensor<2, n, number> Ainv;
  for (std::size_t i = 0; i < n; ++i)
    {
      // i-th column of identity
      dealii::Tensor<1, n, number> e;
      e[i] = number(1.0);

      // Solve L y = e
      dealii::Tensor<1, n, number> y = forward_substitution(L, e);
      // Solve L^T x = y
      dealii::Tensor<1, n, number> x = backward_substitution(L, y);

      // Set i-th column of Ainv
      for (std::size_t j = 0; j < n; ++j)
        {
          Ainv[j][i] = x[j];
        }
    }

  return Ainv;
}
