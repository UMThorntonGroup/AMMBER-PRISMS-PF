#ifndef PARABOLOIDSYSTEM_H
#define PARABOLOIDSYSTEM_H

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <json.hpp>
#include <map>
#include <prismspf/core/field_attributes.h>
#include <prismspf/core/solve_block.h>
#include <string>

using namespace prismspf;

struct Scales
{
  double length_scale, time_scale, energy_scale;

  void
  from_json(const nlohmann::json &j)
  {
    // Parse Dimensions
    length_scale = j.at("dimensions").at("length_scale").get<double>();
    time_scale   = j.at("dimensions").at("time_scale").get<double>();
    energy_scale = j.at("dimensions").at("energy_density_scale").get<double>();
  }
};

/**
 * @brief Class containing the thermodynamic and kinetic parameters needed for the grand
 * potential based model using simplified second-order polynomials
 */
class ParaboloidSystem
{
public:
  /**
   * @brief Parameters for compositions within a phase
   */
  struct PhaseCompInfo
  {
    std::string name;
    double      k_well, c_min, x0;
  };

  /**
   * @brief Parameters for each phase
   */
  struct Phase
  {
    std::string                name;
    double                     mu_int;
    double                     D;
    double                     sigma;
    double                     f_min;
    std::vector<PhaseCompInfo> comps;
  };

  /**
   * @brief Phase parameters
   */
  std::vector<Phase> phases;
  /**
   * @brief Component names at each index
   */
  std::vector<std::string> comp_names;
  /**
   * @brief Phase names at each index
   */
  std::vector<std::string> phase_names;
  /**
   * @brief Name of the component that is not an independent variable
   */
  std::string solution_component;
  /**
   * @brief The index of the phase of each order parameter
   */
  std::vector<uint> order_params;
  /**
   * @brief The atomic or molar volume assumed to be uniform
   */
  double Vm;
  /**
   * @brief The interface width
   */
  double l_int;
  /**
   * @brief If true, the energy density is converted to volumetric energy density from
   * molar/atomic energy
   */
  bool volumetrize;

  /**
   * @brief Constructor
   */
  ParaboloidSystem()
  {}

  /**
   * @brief JSON Constructor
   * @param j JSON object containing the parameters
   */
  ParaboloidSystem(const nlohmann::json &j)
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
    // Parse Vm
    Vm    = j.at("Vm").get<double>();
    l_int = j.at("l_int").get<double>();

    // Check if volumetrization needed
    volumetrize = j.at("convert_fractional_to_volumetric_energy").get<bool>();

    // Parse solution component
    solution_component = j.at("solution_component").get<std::string>();

    // Parse components
    comp_names.clear();
    for (const auto &comp_name : j.at("components"))
      {
        comp_names.push_back(comp_name);
      }

    // Parse phases
    phases.clear();
    phase_names.clear();
    for (const auto &[phase_name, phase_info] : j.at("phases").items())
      {
        Phase phase;
        phase.name = phase_name;
        phase_info.at("mu_int").get_to(phase.mu_int);
        phase_info.at("sigma").get_to(phase.sigma);
        phase_info.at("f_min").get_to(phase.f_min);
        phase_info.at("D").get_to(phase.D);

        // Parse components
        for (const std::string &comp_name : comp_names)
          {
            PhaseCompInfo phaseCompInfo;
            phaseCompInfo.name = comp_name;
            phase_info.at(comp_name).at("c_min").get_to(phaseCompInfo.c_min);
            phase_info.at(comp_name).at("k_well").get_to(phaseCompInfo.k_well);
            phase_info.at(comp_name).at("x0").get_to(phaseCompInfo.x0);
            phase.comps.push_back(phaseCompInfo);
          }
        phases.push_back(phase);
        phase_names.push_back(phase_name);
      }

    // Parse order parameters
    order_params.clear();
    for (const std::string phase_name : j.at("order_parameters"))
      {
        uint phase_index =
          std::find(phase_names.begin(), phase_names.end(), phase_name) - phase_names.begin();
        order_params.push_back(phase_index);
      }

    // Convert to volumetric energy if necessary
    if (volumetrize)
      {
        for (Phase &phase : phases)
          {
            phase.f_min /= Vm;
            for (PhaseCompInfo &comp : phase.comps)
              {
                comp.k_well /= Vm;
              }
          }
      }
  }

  double
  estimate_energy_scale() const
  {
    double minimum_diagonal_curvature = std::numeric_limits<double>::max();
    for (const auto &phase : phases)
      {
        for (const auto &comp : phase.comps)
          {
            minimum_diagonal_curvature = std::min(minimum_diagonal_curvature, comp.k_well);
          }
      }
    return minimum_diagonal_curvature;
  }

  double
  estimate_length_scale() const
  {
    return l_int;
  }

  double
  max_gradient_coefficient() const
  {
    double max_gradient_coefficient = 0.0;
    for (const auto &phase : phases)
      {
        max_gradient_coefficient = std::max(max_gradient_coefficient, phase.D);
        max_gradient_coefficient = std::max(max_gradient_coefficient, phase.mu_int * phase.sigma);
      }
    return max_gradient_coefficient;
  }

  /**
   * @brief Non-dimensionalize the parameters using the provided unit scales.
   */
  void
  nondimensionalize(const Scales &scales)
  {
    const double l0 = scales.length_scale;
    const double t0 = scales.time_scale;
    const double E0 = scales.energy_scale; // energy density
    Vm              = Vm / (l0 * l0 * l0);
    l_int           = l_int / (l0);
    for (Phase &phase : phases)
      {
        phase.mu_int = phase.mu_int / (l0 / (E0 * t0));
        phase.sigma  = phase.sigma / (E0 * l0);
        phase.f_min  = phase.f_min / E0;
        phase.D      = phase.D / ((l0 * l0) / t0);
        for (PhaseCompInfo &comp : phase.comps)
          {
            comp.k_well = comp.k_well / E0;
          }
      }
  }

  /**
   * @brief Print the parameters to the console
   */
  void
  print_parameters()
  {
    // Set column width
    const int col_width  = 15;
    const int line_width = 45;

    // Print header
    std::cout << std::left << std::setw(col_width) << "Name" << std::setw(col_width) << "Dimensionless"
              << std::setw(col_width) << "Dimensional"
              << "\n";
    std::cout << std::string(line_width, '-') << "\n";

    // Print Vm and l_int
    std::cout << std::setw(col_width) << "Vm:" << std::setw(col_width) << Vm << std::setw(col_width) << Vm
              << "\n";
    std::cout << std::setw(col_width) << "l_int:" << std::setw(col_width) << l_int << std::setw(col_width)
              << l_int << "\n";

    // Print component information
    for (const auto &comp_name : comp_names)
      {
        std::cout << std::setw(col_width) << comp_name << "\n";
      }
    std::cout << "\n";

    // Print phase information
    for (const Phase &phase : phases)
      {
        std::cout << std::setw(col_width) << phase.name << "\n";
        std::cout << std::setw(col_width) << "mu_int:" << std::setw(col_width) << phase.mu_int
                  << std::setw(col_width) << phase.mu_int << "\n";
        std::cout << std::setw(col_width) << "D:" << std::setw(col_width) << phase.D << std::setw(col_width)
                  << phase.D << "\n";
        std::cout << std::setw(col_width) << "sigma:" << std::setw(col_width) << phase.sigma
                  << std::setw(col_width) << phase.sigma << "\n";

        for (const PhaseCompInfo &comp : phase.comps)
          {
            std::cout << std::setw(col_width) << comp.name << "\n";
            std::cout << std::setw(col_width) << "k_well:" << std::setw(col_width) << comp.k_well
                      << std::setw(col_width) << comp.k_well << "\n";
            std::cout << std::setw(col_width) << "c_min:" << std::setw(col_width) << comp.c_min
                      << std::setw(col_width) << comp.c_min << "\n";
            std::cout << std::setw(col_width) << "x0:" << std::setw(col_width) << comp.x0 << "\n";
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
    for (unsigned int phase_index = 0; phase_index < phase_names.size(); phase_index++)
      {
        phase_counter[phase_index] = 0;
      }
    for (const unsigned int &phase_index : order_params)
      {
        std::string phase_name = phase_names[phase_index];
        std::string var_name   = phase_name + "_" + std::to_string(phase_counter[phase_index]++);
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
   * @brief Get number of components
   * @return The number of components
   */
  size_t
  num_comps() const
  {
    return comp_names.size();
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
    return mu_base() + num_comps();
  }

  /**
   * @brief  Get base index for order parameters time derivatives
   * @return The base index for order parameters time derivatives
   */
  size_t
  detadt_base() const
  {
    return eta_base() + num_ops();
  }

  /**
   * @brief  Get base index for total concentration
   * @return The base index for total concentration
   */
  size_t
  c_tot_base() const
  {
    return detadt_base() + num_ops();
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
    fields.reserve(2 * comp_names.size() + 2 * op_names.size());
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

    for (uint var_index = mu_base(); var_index < num_comps(); var_index++)
      {
        mu_and_eta.field_indices.insert(var_index);
      }
    for (uint var_index = eta_base(); var_index < num_ops(); var_index++)
      {
        mu_and_eta.field_indices.insert(var_index);
      }
    for (uint var_index = detadt_base(); var_index < detadt_base() + num_ops(); var_index++)
      {
        detadt.field_indices.insert(var_index);
      }
    for (uint var_index = c_tot_base(); var_index < c_tot_base() + num_comps(); var_index++)
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
};

#endif