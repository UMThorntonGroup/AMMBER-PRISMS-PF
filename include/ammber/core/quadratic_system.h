#pragma once

#include <ammber/core/scales.h>

#include <prismspf/core/field_attributes.h>
#include <prismspf/core/solve_block.h>

#include <deal.II/base/tensor.h>
#include <deal.II/lac/full_matrix.h>

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

  std::vector<std::string>   comp_names;
  std::vector<Phase>         phases;
  std::vector<std::string>   phase_names;
  std::string                solution_component;
  std::vector<unsigned int>  order_params;
  dealii::FullMatrix<double> sigma;
  double                     Vm    = 0.0;
  double                     l_int = 0.0;

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
    if (comp_names.size() != num_comps)
      {
        throw std::invalid_argument("Quadratic system component count does not match its template");
      }

    phases.clear();
    phase_names.clear();
    for (const auto &phase_json : j.at("phases"))
      {
        Phase phase;
        phase.name   = phase_json.at("name").get<std::string>();
        phase.mu_int = phase_json.at("mu_int").get<double>();
        phase.C      = phase_json.at("C").get<double>();
        phase.D      = parse_component_matrix(phase_json.at("D"));
        phase.A      = parse_component_matrix(phase_json.at("A"));
        phase.B      = parse_vector(phase_json.at("B"));
        phase.c_ref  = parse_vector(phase_json.at("c_ref"));
        phase.c0     = parse_vector(phase_json.at("c0"));

        if (std::find(phase_names.begin(), phase_names.end(), phase.name) != phase_names.end())
          {
            throw std::invalid_argument("Duplicate phase name in quadratic system: " + phase.name);
          }
        phase_names.push_back(phase.name);
        phases.push_back(std::move(phase));
      }

    sigma = parse_sigma_matrix(j.at("sigma"), phases.size());

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
        const auto print_matrix = [&](const char *name, const CompMatrix<double> &matrix)
        {
          std::cout << name << "\n";
          for (unsigned int i = 0; i < num_comps; ++i)
            {
              for (unsigned int j = 0; j < num_comps; ++j)
                {
                  std::cout << std::setw(col_width) << matrix[i][j];
                }
              std::cout << "\n";
            }
        };
        print_matrix("D:", phase.D);
        print_matrix("A:", phase.A);
        std::cout << "B:";
        for (unsigned int i = 0; i < num_comps; ++i)
          {
            std::cout << std::setw(col_width) << phase.B[i];
          }
        std::cout << "\nc_ref:";
        for (unsigned int i = 0; i < num_comps; ++i)
          {
            std::cout << std::setw(col_width) << phase.c_ref[i];
          }
        std::cout << "\nc0:";
        for (unsigned int i = 0; i < num_comps; ++i)
          {
            std::cout << std::setw(col_width) << phase.c0[i];
          }
        std::cout << "\nC:" << std::setw(col_width) << phase.C << "\n\n";
      }

    std::cout << "sigma:\n";
    for (std::size_t i = 0; i < sigma.size()[0]; ++i)
      {
        for (std::size_t j = 0; j < sigma.size()[1]; ++j)
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
        phase.mu0 = phase.A * (phase.c0 - phase.c_ref) + phase.B;

        dealii::FullMatrix<double> helper(num_comps);
        for (unsigned int i = 0; i < num_comps; ++i)
          {
            for (unsigned int j = 0; j < num_comps; ++j)
              {
                helper(i, j) = phase.A[i][j];
              }
          }
        dealii::FullMatrix<double> helper_inv(num_comps);
        helper_inv.invert(helper);
        for (unsigned int i = 0; i < num_comps; ++i)
          {
            for (unsigned int j = 0; j < num_comps; ++j)
              {
                phase.suscept[i][j] = helper_inv(i, j);
              }
          }
      }
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
    return num_comps;
    // return mu_base() + num_comps;
  }

  /**
   * @brief  Get base index for order parameters time derivatives
   * @return The base index for order parameters time derivatives
   */
  size_t
  detadt_base() const
  {
    return num_comps + num_ops();
    // return eta_base() + num_ops();
  }

  /**
   * @brief  Get base index for total concentration
   * @return The base index for total concentration
   */
  size_t
  c_tot_base() const
  {
    return num_comps + 2 * num_ops();
    // return detadt_base() + num_ops();
  }

  /**
   * @brief  Get base index for phase fraction h
   * @return The base index for phase fraction h
   */
  size_t
  h_base() const
  {
    return 2 * num_comps + 2 * num_ops();
    // return c_tot_base() + num_comps;
  }

  /**
   * @brief  Get base index for grand potential density omega
   * @return The base index for grand potential density omega
   */
  size_t
  omega_base() const
  {
    return 2 * num_comps + 2 * num_ops() + num_phases();
    // return h_base() + num_phases();
  }

  static constexpr unsigned int explicit_block_id = 0;
  static constexpr unsigned int detadt_block_id   = 1;
  static constexpr unsigned int pp_block_id       = 2;

  /**
   * @brief Declare the fields needed for the PDE in PRISMS-PF
   * @param loader Pointer to the attribute loader (equations.cc)
   */
  std::vector<FieldAttributes>
  load_fields()
  {
    // Get names for order parameter fields
    std::vector<std::string> op_names = get_order_parameter_names();

    std::vector<FieldAttributes> fields;
    fields.reserve(2 * num_comps + 2 * num_ops() + 2 * num_phases());
    for (const auto &comp_name : comp_names)
      {
        fields.emplace_back("mu_" + comp_name);
      }
    for (const auto &op_name : op_names)
      {
        fields.emplace_back("eta_" + op_name);
      }
    for (const auto &op_name : op_names)
      {
        fields.emplace_back("detadt_" + op_name);
      }
    for (const auto &comp_name : comp_names)
      {
        fields.emplace_back("c_" + comp_name);
      }
    for (const auto &phase_name : phase_names)
      {
        fields.emplace_back("h_" + phase_name);
      }
    for (const auto &phase_name : phase_names)
      {
        fields.emplace_back("omega_" + phase_name);
      }
    return fields;
  }

  /**
   * @brief Declare the fields needed for the PDE in PRISMS-PF
   * @param loader Pointer to the attribute loader (equations.cc)
   */
  std::vector<SolveBlock>
  load_blocks()
  {
    static const Dependency old_1_val_and_grad(EvalFlags::nothing,
                                               EvalFlags::nothing,
                                               {EvalFlags::values | EvalFlags::gradients});
    static const Dependency current_val_and_grad(EvalFlags::values | EvalFlags::gradients);
    static const Dependency current_val(EvalFlags::values);

    SolveBlock mu_and_eta;
    SolveBlock detadt;
    SolveBlock pp;

    for (uint var_index = mu_base(); var_index < mu_base() + num_comps; var_index++)
      {
        mu_and_eta.field_indices.insert(var_index);
      }
    for (uint var_index = eta_base(); var_index < eta_base() + num_ops(); var_index++)
      {
        mu_and_eta.field_indices.insert(var_index);
      }
    for (uint var_index = detadt_base(); var_index < detadt_base() + num_ops(); var_index++)
      {
        detadt.field_indices.insert(var_index);
      }
    for (uint var_index = c_tot_base(); var_index < c_tot_base() + num_comps; var_index++)
      {
        pp.field_indices.insert(var_index);
      }
    for (uint var_index = h_base(); var_index < h_base() + num_phases(); var_index++)
      {
        pp.field_indices.insert(var_index);
      }
    for (uint var_index = omega_base(); var_index < omega_base() + num_phases(); var_index++)
      {
        pp.field_indices.insert(var_index);
      }

    mu_and_eta.id           = explicit_block_id;
    mu_and_eta.solve_type   = Explicit;
    mu_and_eta.solve_timing = Primary;
    for (uint field_index : mu_and_eta.field_indices)
      {
        mu_and_eta.dependencies_rhs[field_index] = old_1_val_and_grad;
      }
    for (uint field_index : detadt.field_indices)
      {
        mu_and_eta.dependencies_rhs[field_index] = old_1_val_and_grad;
      }

    detadt.id           = detadt_block_id;
    detadt.solve_type   = Explicit;
    detadt.solve_timing = Secondary;
    for (uint field_index : mu_and_eta.field_indices)
      {
        detadt.dependencies_rhs[field_index] = current_val_and_grad;
      }

    pp.id           = pp_block_id;
    pp.solve_type   = Explicit;
    pp.solve_timing = PostProcess;
    for (uint field_index : mu_and_eta.field_indices)
      {
        pp.dependencies_rhs[field_index] = current_val;
      }

    return std::vector<SolveBlock>({mu_and_eta, detadt, pp});
  }

  /**
   * @brief Function to set the initial conditions with the proper tanh profile given a
   * level-set function
   */
  template <typename number>
  [[nodiscard]] auto
  interface(const number &sdf) const
  {
    using std::tanh;
    return 0.5 * (1.0 + tanh(2.0 * sdf / l_int));
  }

private:
  static CompVector<double>
  parse_vector(const nlohmann::json &j)
  {
    if (!j.is_array() || j.size() != num_comps)
      {
        throw std::invalid_argument("Quadratic phase vector has an unexpected size");
      }

    CompVector<double> result;
    for (unsigned int i = 0; i < num_comps; ++i)
      {
        result[i] = j.at(i).get<double>();
      }
    return result;
  }

  static CompMatrix<double>
  parse_component_matrix(const nlohmann::json &j)
  {
    if (!j.is_array() || j.size() != num_comps)
      {
        throw std::invalid_argument("Quadratic system matrix has an unexpected row count");
      }

    CompMatrix<double> result;
    for (unsigned int i = 0; i < num_comps; ++i)
      {
        if (!j.at(i).is_array() || j.at(i).size() != num_comps)
          {
            throw std::invalid_argument("Quadratic system matrix has an unexpected column count");
          }
        for (unsigned int k = 0; k < num_comps; ++k)
          {
            result[i][k] = j.at(i).at(k).get<double>();
          }
      }
    return result;
  }

  static dealii::FullMatrix<double>
  parse_sigma_matrix(const nlohmann::json &j, const std::size_t size)
  {
    if (!j.is_array() || j.size() != size)
      {
        throw std::invalid_argument("Quadratic system matrix has an unexpected row count");
      }

    dealii::FullMatrix<double> result(size, size);
    for (std::size_t i = 0; i < size; ++i)
      {
        if (!j.at(i).is_array() || j.at(i).size() != size)
          {
            throw std::invalid_argument("Quadratic system matrix has an unexpected column count");
          }
        for (std::size_t k = 0; k < size; ++k)
          {
            result(i, k) = j.at(i).at(k).get<double>();
          }
      }
    return result;
  }
};