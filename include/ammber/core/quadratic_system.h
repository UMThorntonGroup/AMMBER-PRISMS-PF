#pragma once

#include <ammber/core/scales.h>

#include <prismspf/core/field_attributes.h>
#include <prismspf/core/solve_block.h>

#include <deal.II/lac/full_matrix.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace prismspf;

/**
 * @brief Class containing the thermodynamic and kinetic parameters needed for the grand
 * potential based model using second-order polynomials
 */
template <unsigned num_comps>
class QuadraticSystem
{
public:
  template <typename T>
  using CompVector = dealii::Tensor<1, num_comps, T>;
  template <typename T>
  using CompMatrix = dealii::Tensor<2, num_comps, T>;
  template <typename T>
  using CompSymmet = dealii::SymmetricTensor<2, num_comps, T>;

  /**
   * @brief Parameters for each phase
   */
  struct Phase
  {
    std::string        name;
    double             mu_int = 0.0;
    CompMatrix<double> D, A, suscept;
    CompVector<double> B, c_ref, c0, mu0;
    double             C = 0.0;
  };

  std::array<std::string, num_comps> comp_names;
  std::vector<Phase>                 phases;
  std::vector<std::string>           phase_names;
  std::string                        solution_component;
  std::vector<unsigned int>          order_params;
  dealii::FullMatrix<double>         sigma;
  double                             Vm    = 0.0;
  double                             l_int = 0.0;

  /**
   * @brief Constructor
   */
  QuadraticSystem() = default;

  /**
   * @brief JSON Constructor
   * @param j JSON object containing the parameters
   */
  explicit QuadraticSystem(const nlohmann::json &j)
  {
    from_json(j);
  }

  /**
   * @brief Load the parameters from a JSON object
   * @param j JSON object containing the parameters
   */
  void
  from_json(const nlohmann::json &j)
  {
    Vm                     = j.at("Vm").get<double>();
    l_int                  = j.at("l_int").get<double>();
    const bool volumetrize = j.at("convert_fractional_to_volumetric_energy").get<bool>();
    solution_component     = j.at("solution_component").get<std::string>();

    comp_names = j.at("components").get<std::vector<std::string>>();
    num_comps  = static_cast<unsigned int>(comp_names.size());

    phases.clear();
    phase_names.clear();
    for (const auto &phase_json : j.at("phases"))
      {
        Phase phase;
        phase.name   = phase_json.at("name").get<std::string>();
        phase.mu_int = phase_json.at("mu_int").get<double>();
        phase.C      = phase_json.at("C").get<double>();
        phase.D      = parse_matrix(phase_json.at("D"), num_comps, num_comps);
        phase.A      = parse_matrix(phase_json.at("A"), num_comps, num_comps);
        phase.B      = parse_vector(phase_json.at("B"), num_comps);
        phase.c_ref  = parse_vector(phase_json.at("c_ref"), num_comps);
        phase.c0     = parse_vector(phase_json.at("c0"), num_comps);

        if (std::find(phase_names.begin(), phase_names.end(), phase.name) != phase_names.end())
          {
            throw std::invalid_argument("Duplicate phase name in quadratic system: " + phase.name);
          }
        phase_names.push_back(phase.name);
        phases.push_back(std::move(phase));
      }

    sigma = parse_matrix(j.at("sigma"), phases.size(), phases.size());

    order_params.clear();
    for (const auto &phase_name_json : j.at("order_parameters"))
      {
        const std::string phase_name = phase_name_json.get<std::string>();
        const auto        phase_iter = std::find(phase_names.begin(), phase_names.end(), phase_name);
        if (phase_iter == phase_names.end())
          {
            throw std::invalid_argument("Unknown order parameter phase: " + phase_name);
          }
        order_params.push_back(static_cast<uint>(phase_iter - phase_names.begin()));
      }

    if (volumetrize)
      {
        if (Vm <= 0.0)
          {
            throw std::invalid_argument(
              "Vm must be positive to convert fractional energy to volumetric energy");
          }
        for (Phase &phase : phases)
          {
            phase.A /= Vm;
            phase.B /= Vm;
            phase.C /= Vm;
          }
      }

    update_derived_parameters();
  }

  /**
   * @brief Non-dimensionalize the parameters using the provided unit scales
   */
  void
  nondimensionalize(const Scales &scales)
  {
    if (scales.length_scale <= 0.0 || scales.time_scale <= 0.0 || scales.energy_scale <= 0.0)
      {
        throw std::invalid_argument(
          "Quadratic system dimensions must all be positive before nondimensionalizing");
      }

    const double l0 = scales.length_scale;
    const double t0 = scales.time_scale;
    const double E0 = scales.energy_scale;
    Vm /= (l0 * l0 * l0);
    l_int /= l0;
    sigma /= (E0 * l0);
    for (Phase &phase : phases)
      {
        phase.mu_int /= (l0 / (E0 * t0));
        phase.D /= ((l0 * l0) / t0);
        phase.A /= E0;
        phase.B /= E0;
        phase.C /= E0;
      }
    update_derived_parameters();
  }

  /**
   * @brief Print the parameters to the console
   */
  void
  print_parameters() const
  {
    const int col_width = 15;
    std::cout << std::left << std::setw(col_width) << "Name" << std::setw(col_width) << "Value"
              << "\n";
    std::cout << std::string(2 * col_width, '-') << "\n";
    std::cout << std::setw(col_width) << "Vm:" << std::setw(col_width) << Vm << "\n";
    std::cout << std::setw(col_width) << "l_int:" << std::setw(col_width) << l_int << "\n\n";

    for (const Phase &phase : phases)
      {
        std::cout << std::setw(col_width) << phase.name << "\n";
        std::cout << std::setw(col_width) << "mu_int:" << std::setw(col_width) << phase.mu_int << "\n";
        const auto print_matrix = [&](const char *name, const boost_matrix<double> &matrix)
          {
            std::cout << name << "\n";
            for (std::size_t i = 0; i < matrix.size1(); ++i)
              {
                for (std::size_t j = 0; j < matrix.size2(); ++j)
                  {
                    std::cout << std::setw(col_width) << matrix(i, j);
                  }
                std::cout << "\n";
              }
          };
        print_matrix("D:", phase.D);
        print_matrix("A:", phase.A);
        std::cout << "B:";
        for (std::size_t i = 0; i < phase.B.size(); ++i)
          {
            std::cout << std::setw(col_width) << phase.B[i];
          }
        std::cout << "\nc_ref:";
        for (std::size_t i = 0; i < phase.c_ref.size(); ++i)
          {
            std::cout << std::setw(col_width) << phase.c_ref[i];
          }
        std::cout << "\nc0:";
        for (std::size_t i = 0; i < phase.c0.size(); ++i)
          {
            std::cout << std::setw(col_width) << phase.c0[i];
          }
        std::cout << "\nC:" << std::setw(col_width) << phase.C << "\n\n";
      }

    std::cout << "sigma:\n";
    for (std::size_t i = 0; i < sigma.size1(); ++i)
      {
        for (std::size_t j = 0; j < sigma.size2(); ++j)
          {
            std::cout << std::setw(col_width) << sigma(i, j);
          }
        std::cout << "\n";
      }
  }

  /**
   * @brief Get names of order parameters
   * @return A vector of order parameter names
   */
  std::vector<std::string>
  get_order_parameter_names() const
  {
    std::vector<std::string> op_names;
    std::map<uint, uint>     phase_counter;
    for (uint phase_index : order_params)
      {
        const std::string var_name =
          phase_names.at(phase_index) + "_" + std::to_string(phase_counter[phase_index]++);
        op_names.push_back(var_name);
      }
    return op_names;
  }

  /**
   * @brief Get names of components
   * @return A vector of component names
   */
  const std::vector<std::string> &
  get_component_names() const
  {
    return comp_names;
  }

  /**
   * @brief Declare the fields needed for the PDE in PRISMS-PF
   * @param loader Pointer to the attribute loader (equations.cc)
   */
  void
  load_variables()
  {}

  void
  update_derived_parameters()
  {
    for (Phase &phase : phases)
      {
        phase.mu0 = prod(phase.A, phase.c0 - phase.c_ref) + phase.B;

        dealii::FullMatrix<double> helper(num_comps);
        for (std::size_t i = 0; i < num_comps; ++i)
          {
            for (std::size_t j = 0; j < num_comps; ++j)
              {
                helper(i, j) = phase.A[i][j];
              }
          }
        dealii::FullMatrix<double> helper_inv(num_comps);
        helper_inv.invert(helper);
        phase.suscept = boost_matrix<double>(num_comps, num_comps);
        for (std::size_t i = 0; i < num_comps; ++i)
          {
            for (std::size_t j = 0; j < num_comps; ++j)
              {
                phase.suscept[i][j] = helper_inv(i, j);
              }
          }
      }
  }

  /**
   * @brief Get number of components
   * @return The number of components
   */
  size_t
  num_comps() const
  {
    return comp_names.size();
  }

  /**
   * @brief  Get number of phases
   * @return The number of phases
   */
  size_t
  num_phases() const
  {
    return phases.size();
  }

  /**
   * @brief Get number of order parameters
   * @return The number of order parameters
   */
  size_t
  num_ops() const
  {
    return order_params.size();
  }

  /**
   * @brief Get base index for diffusion potentials
   * @return The base index for diffusion potentials
   */
  size_t
  mu_base() const
  {
    return 0;
  }

  /**
   * @brief Get base index for order parameters
   * @return The base index for order parameters
   */
  size_t
  eta_base() const
  {
    return num_comps();
    // return mu_base() + num_comps();
  }

  /**
   * @brief  Get base index for phase fraction h
   * @return The base index for phase fraction h
   */
  size_t
  h_base() const
  {
    return 2 * num_comps() + 2 * num_ops();
    // return c_tot_base() + num_comps();
  }

  /**
   * @brief  Get base index for order parameters time derivatives
   * @return The base index for order parameters time derivatives
   */
  size_t
  detadt_base() const
  {
    return num_comps() + num_ops();
    // return eta_base() + num_ops();
  }

  /**
   * @brief  Get base index for total concentration
   * @return The base index for total concentration
   */
  size_t
  c_tot_base() const
  {
    return num_comps() + 2 * num_ops();
    // return detadt_base() + num_ops();
  }

  /**
   * @brief  Get base index for grand potential density omega
   * @return The base index for grand potential density omega
   */
  size_t
  omega_base() const
  {
    return 2 * num_comps() + 2 * num_ops() + num_phases();
    // return h_base() + num_phases();
  }

private:
  static boost_vector<double>
  parse_vector(const nlohmann::json &j, std::size_t size)
  {
    if (!j.is_array() || j.size() != size)
      {
        throw std::invalid_argument("Quadratic phase vector has an unexpected size");
      }

    boost_vector<double> result(size);
    for (std::size_t i = 0; i < size; ++i)
      {
        result[i] = j.at(i).get<double>();
      }
    return result;
  }

  static boost_matrix<double>
  parse_matrix(const nlohmann::json &j, std::size_t rows, std::size_t columns)
  {
    if (!j.is_array() || j.size() != rows)
      {
        throw std::invalid_argument("Quadratic system matrix has an unexpected row count");
      }

    boost_matrix<double> result(rows, columns);
    for (std::size_t i = 0; i < rows; ++i)
      {
        if (!j.at(i).is_array() || j.at(i).size() != columns)
          {
            throw std::invalid_argument("Quadratic system matrix has an unexpected column count");
          }
        for (std::size_t k = 0; k < columns; ++k)
          {
            result(i, k) = j.at(i).at(k).get<double>();
          }
      }
    return result;
  }
};