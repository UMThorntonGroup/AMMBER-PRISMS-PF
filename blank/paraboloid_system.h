#pragma once

#include <iomanip>
#include <iostream>
#include <json.hpp>
#include <map>
#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/user_inputs/user_input_parameters.h>
#include <string>

struct Scales
{
  double length_scale;
  double energy_scale;
  double time_scale;

  Scales()
    : length_scale(0.0)
    , energy_scale(0.0)
    , time_scale(0.0)
  {}

  Scales(nlohmann::json j)
    : length_scale(j.at("length_scale").get<double>())
    , energy_scale(j.at("energy_density_scale").get<double>())
    , time_scale(j.at("time_scale").get<double>())
  {}
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
   * @brief Whether to convert fractional free energies to volumetric free energies
   */
  bool volumetrize;
  /**
   * @brief The interface width
   */
  double l_int;

  /**
   * @brief Constructor
   */
  ParaboloidSystem()
  {}

  /**
   * @brief JSON Constructor
   * @param TCSystem JSON object containing the parameters
   */
  ParaboloidSystem(const nlohmann::json &TCSystem)
  {
    from_json(TCSystem);
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
        uint phase_index = std::find(phase_names.begin(), phase_names.end(), phase_name) -
                           phase_names.begin();
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

  /**
   * @brief Automatically select scales if they are set to zero
   */
  template <unsigned int dim>
  Scales
  select_scales(const Scales                           &input_scales,
                const prisms::UserInputParameters<dim> &user_inputs) const
  {
    Scales                                    output_scales = input_scales;
    const prisms::SpatialDiscretization<dim> &spc_disc =
      user_inputs.get_spatial_discretization();
    const prisms::TemporalDiscretization &temp_disc =
      user_inputs.get_temporal_discretization();
    const prisms::UserConstants<dim> &user_constants = user_inputs.get_user_constants();
    // If energy_scale is zero, set it to the lowest free energy curvature
    // ========================================================================
    if (input_scales.energy_scale == 0.0)
      {
        double minimum_diagonal_curvature = std::numeric_limits<double>::max();
        for (const auto &phase : phases)
          {
            for (const auto &comp : phase.comps)
              {
                minimum_diagonal_curvature =
                  std::min(minimum_diagonal_curvature, comp.k_well);
              }
          }
        output_scales.energy_scale = minimum_diagonal_curvature;
        prisms::ConditionalOStreams::pout_base()
          << "Setting energy scale to " << output_scales.energy_scale
          << " based on the lowest free energy curvature diagonal.\n";
      }
    // ========================================================================
    // If l_int is zero, guess it based on thermodynamics
    // ========================================================================
    // if (_l_int == 0.0)
    //   {
    //     _l_int = phases[0]._f_min;
    //   }
    // ========================================================================

    // If length_scale is zero, set it based on the interface width and domain
    // discretization
    // ========================================================================
    double min_dx = std::numeric_limits<double>::max();
    double points_in_interface =
      user_constants.get_model_constant_double("num_points_in_interface");
    for (unsigned int i = 0; i < dim; i++)
      {
        min_dx =
          std::min(min_dx,
                   double(spc_disc.get_subdivisions()[i]) * spc_disc.get_size()[i] /
                     std::pow(2.0,
                              std::max(spc_disc.get_global_refinement(),
                                       spc_disc.get_max_refinement())));
      }
    if (input_scales.length_scale == 0.0)
      {
        output_scales.length_scale = l_int / (min_dx * points_in_interface);
        prisms::ConditionalOStreams::pout_base()
          << "Setting length scale to " << output_scales.length_scale
          << " based on the interface width and domain discretization using "
          << points_in_interface << " points in the interface.\n";
      }
    // ========================================================================

    const double time_scale_factor = user_constants.get_model_constant_double(
      "time_scale_stability_factor"); // Design factor for stability
    constexpr double theoretical_max_gradient_factor = 0.25;
    double           max_gradient_factor             = 0.0;
    std::string      stability_limiter               = "diffusion";
    std::string      limiting_phase;
    const double     gradient_prefactor =
      temp_disc.get_timestep() * (spc_disc.get_degree() * spc_disc.get_degree()) /
      (min_dx * min_dx * output_scales.length_scale * output_scales.length_scale);
    for (const auto &phase : phases)
      {
        const double diffusion_gradient_factor = phase.D * gradient_prefactor;
        // mu*sigma = L*kappa
        const double order_parameter_gradient_factor =
          (phase.mu_int * phase.sigma) * gradient_prefactor;

        if (diffusion_gradient_factor > max_gradient_factor)
          {
            max_gradient_factor = diffusion_gradient_factor;
            stability_limiter   = "diffusion";
            limiting_phase      = phase.name;
          }
        if (order_parameter_gradient_factor > max_gradient_factor)
          {
            max_gradient_factor = order_parameter_gradient_factor;
            stability_limiter   = "order parameter evolution";
            limiting_phase      = phase.name;
          }
      }
    // If time_scale is zero, set it based on the stability limit
    // ========================================================================
    const double reccommended_time_scale =
      time_scale_factor * theoretical_max_gradient_factor / max_gradient_factor;
    prisms::ConditionalOStreams::pout_base()
      << "The numerical stability for this set of parameters is limited by "
      << stability_limiter << " in phase " << limiting_phase << ".\n";
    if (input_scales.time_scale == 0.0)
      {
        output_scales.time_scale = reccommended_time_scale;
        prisms::ConditionalOStreams::pout_base()
          << "\nSetting time scale to " << output_scales.time_scale
          << " based on the stability limit using a design factor of "
          << time_scale_factor << ".\n";
      }
    else if (input_scales.time_scale > reccommended_time_scale)
      { // If time_scale is not zero, ensure it is not larger than the stability limit
        prisms::ConditionalOStreams::pout_base()
          << "\nWarning: Provided time scale is larger than the stability "
             "limit using a design factor of "
          << time_scale_factor
          << ".\n"
             "We recommend using a time scale of "
          << reccommended_time_scale << " or a time step of "
          << temp_disc.get_timestep() * reccommended_time_scale / input_scales.time_scale
          << " to ensure stability.\n\n";
      }
    return output_scales;
    // ========================================================================
  }

  /**
   * @brief Print the parameters to the console
   */
  void
  print_parameters() const
  {
    // Set column width
    const int col_width  = 15;
    const int line_width = 45;

    // Print header
    prisms::ConditionalOStreams::pout_base()
      << std::left << std::setw(col_width) << "Name" << std::setw(col_width) << "Value"
      << "\n";
    prisms::ConditionalOStreams::pout_base() << std::string(line_width, '-') << "\n";

    // Print Vm and l_int
    prisms::ConditionalOStreams::pout_base()
      << std::setw(col_width) << "Vm:" << std::setw(col_width) << Vm << "\n";
    prisms::ConditionalOStreams::pout_base()
      << std::setw(col_width) << "l_int:" << std::setw(col_width) << l_int << "\n";

    // Print component information
    for (const auto &comp_name : comp_names)
      {
        prisms::ConditionalOStreams::pout_base()
          << std::setw(col_width) << comp_name << "\n";
      }
    prisms::ConditionalOStreams::pout_base() << "\n";

    // Print phase information
    for (const Phase &phase : phases)
      {
        prisms::ConditionalOStreams::pout_base()
          << std::setw(col_width) << phase.name << "\n";
        prisms::ConditionalOStreams::pout_base()
          << std::setw(col_width) << "mu_int:" << std::setw(col_width) << phase.mu_int
          << "\n";
        prisms::ConditionalOStreams::pout_base()
          << std::setw(col_width) << "D:" << std::setw(col_width) << phase.D << "\n";
        prisms::ConditionalOStreams::pout_base()
          << std::setw(col_width) << "sigma:" << std::setw(col_width) << phase.sigma
          << "\n";

        for (const PhaseCompInfo &comp : phase.comps)
          {
            prisms::ConditionalOStreams::pout_base()
              << std::setw(col_width) << comp.name << "\n";
            prisms::ConditionalOStreams::pout_base()
              << std::setw(col_width) << "k_well:" << std::setw(col_width) << comp.k_well
              << "\n";
            prisms::ConditionalOStreams::pout_base()
              << std::setw(col_width) << "c_min:" << std::setw(col_width) << comp.c_min
              << "\n";
            prisms::ConditionalOStreams::pout_base()
              << std::setw(col_width) << "x0:" << std::setw(col_width) << comp.x0 << "\n";
          }
        prisms::ConditionalOStreams::pout_base() << "\n";
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
    for (const auto &phase_index : order_params)
      {
        if (phase_counter.find(phase_index) == phase_counter.end())
          {
            phase_counter[phase_index] = 0;
          }
        std::string phase_name = phase_names[phase_index];
        std::string var_name =
          phase_name + "_" + std::to_string(phase_counter[phase_index]);
        op_names.push_back(var_name);
        phase_counter[phase_index]++;
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
   * @brief
   * @return
   */
  template <typename number>
  inline const auto
  interface(number x) const
  {
    using std::tanh;
    return 0.5 * (1.0 + tanh(2.0 * x / l_int));
  }
};

/**
 * @brief Non-dimensionalize the parameters using the provided unit scales.
 */
inline ParaboloidSystem
nondimensionalize(const ParaboloidSystem &sys, const Scales &scales)
{
  ParaboloidSystem nondim_system = sys;
  const double    &l0            = scales.length_scale;
  const double    &t0            = scales.time_scale;
  const double    &E0            = scales.energy_scale; // energy density
  nondim_system.Vm /= (l0 * l0 * l0);
  nondim_system.l_int /= (l0);
  for (ParaboloidSystem::Phase &phase : nondim_system.phases)
    {
      phase.mu_int /= (l0 / (E0 * t0));
      phase.sigma /= (E0 * l0);
      phase.f_min /= E0;
      phase.D /= ((l0 * l0) / t0);
      for (ParaboloidSystem::PhaseCompInfo &comp : phase.comps)
        {
          comp.k_well /= E0;
        }
    }
  return nondim_system;
}