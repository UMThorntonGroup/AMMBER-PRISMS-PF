#include "system_equations.h"

#include <string>

template <int dim, int degree>
class customPDE : public MatrixFreePDE<dim, degree>
{
public:
  // Constructor
  customPDE(UserInputParameters<dim> _userInputs)
    : MatrixFreePDE<dim, degree>(_userInputs)
    , userInputs(_userInputs)
  {
    // Load the model parameters
    std::ifstream ifs("system.json");
    if (!ifs.is_open())
      {
        throw std::runtime_error("Could not open system.json");
      }
    nlohmann::json model_parameters;
    ifs >> model_parameters;
    isoSys.from_json(model_parameters, &userInputs);
    isoSys.print_parameters();
    print_initial_energies();
    print_interface_properties();
  }

  /**
   * @brief Function to set the initial conditions with the proper tanh profile given a
   * level-set function
   */
  [[nodiscard]] double
  interface(double x) const
  {
    return 0.5 * (1.0 + std::tanh(2.0 * x / isoSys.l_int));
  }

  /**
   * @brief Prints the grand potential densities (nondimensionalized) of each phase at
   * its initial conditions
   */
  void
  print_initial_energies()
  {
    std::cout << "Initial omega free energies:\n";
    for (uint phase_index = 0; phase_index < isoSys.phases.size(); phase_index++)
      {
        const ParaboloidSystem::Phase &phase = isoSys.phases.at(phase_index);
        SystemContainer<dim, degree>   sys_for_print(isoSys, userInputs);
        sys_for_print.op_data.push_back({
          phase_index,
          {
            {constV(1.), {}}, // eta
            {constV(0.), {}}, // detadt
            constV(0.), // detadt_field
            {}                // dhdeta
          }
        });

        for (uint comp_index = 0; comp_index < isoSys.comp_names.size(); comp_index++)
          {
            const ParaboloidSystem::PhaseCompInfo &comp_info = phase.comps.at(comp_index);
            double mu0 = comp_info.k_well * (comp_info.x0 - comp_info.c_min);
            sys_for_print.comp_data[comp_index].mu.val = constV(mu0);
          }
        sys_for_print.calculate_omega_phase();

        std::cout << phase.name << ":\n"
                  << "Omega:\t" << sys_for_print.phase_data[phase_index].omega.val[0]
                  << "\n";
        initial_omega[phase.name] = sys_for_print.phase_data[phase_index].omega.val[0];
        for (uint comp_index = 0; comp_index < isoSys.comp_names.size(); comp_index++)
          {
            const ParaboloidSystem::PhaseCompInfo &comp = phase.comps.at(comp_index);
            std::cout << "mu_" << comp.name << ":\t"
                      << sys_for_print.comp_data[comp_index].mu.val[0] << "\n";
          }
        std::cout << "\n";
      }
  }

  /**
   * @brief Function to print the properties of the interface between two phases at
   * initial conditions
   */
  void
  print_interface_properties()
  {
    for (const auto &alpha : isoSys.phases)
      {
        for (const auto &beta : isoSys.phases)
          {
            std::cout << "Properties of the interface between " << alpha.name << " and "
                      << beta.name << ":\n";
            double delta_g = initial_omega[alpha.name] - initial_omega[beta.name];
            double sigma   = 0.5 * (alpha.sigma + beta.sigma);
            double D       = 0.5 * (alpha.D * beta.D) / (alpha.D + beta.D);
            double mu_int =
              0.5 * (alpha.mu_int * beta.mu_int) / (alpha.mu_int + beta.mu_int);
            std::cout << "The value of the dimensionless number delta_g/(sigma/l_int) is "
                      << delta_g * isoSys.l_int / sigma << "\n";
            std::cout
              << "The value of the dimensionless number delta_g*mu_int*l_int/D is "
              << delta_g * mu_int * isoSys.l_int / D << "\n";
            std::cout << "The value of the dimensionless number sigma*mu_int/D is "
                      << sigma * mu_int / D << "\n\n";
          }
      }
  }

  // ================================================================
  // Members specific to this subclass
  // ================================================================
  /**
   * @brief JSON containing the model parameters
   */
  nlohmann::json model_parameters;
  /**
   * @brief Object containing the thermodynamic and kinetic parameters
   */
  ParaboloidSystem isoSys;

  /**
   * @brief Map of the initial grand potential densities for each phase (used for
   * printing)
   */
  std::map<std::string, double> initial_omega;

  double r0 = userInputs.get_model_constant_double("r0");
  // ================================================================
};

// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/pde_operator_base.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class CustomPDE : public PDEOperatorBase<dim, degree, number>
{
public:
  using ScalarValue = dealii::VectorizedArray<number>;
  using ScalarGrad  = dealii::Tensor<1, dim, ScalarValue>;
  using ScalarHess  = dealii::Tensor<2, dim, ScalarValue>;
  using VectorValue = dealii::Tensor<1, dim, ScalarValue>;
  using VectorGrad  = dealii::Tensor<2, dim, ScalarValue>;
  using VectorHess  = dealii::Tensor<3, dim, ScalarValue>;
  using PDEOperatorBase<dim, degree, number>::get_user_inputs;
  using PDEOperatorBase<dim, degree, number>::get_pf_tools;

  ParaboloidSystem sys;

  /**
   * @brief Constructor.
   */
  CustomPDE(const UserInputParameters<dim> &_user_inputs, PhaseFieldTools<dim> &_pf_tools)
    : PDEOperatorBase<dim, degree, number>(_user_inputs, _pf_tools)
  {}

private:
  void
  set_initial_condition([[maybe_unused]] const unsigned int       &index,
                        [[maybe_unused]] const unsigned int       &component,
                        [[maybe_unused]] const dealii::Point<dim> &point,
                        [[maybe_unused]] number                   &scalar_value,
                        [[maybe_unused]] number &vector_component_value) const override
  {
    const dealii::Tensor<1, dim> &mesh_size =
      get_user_inputs().spatial_discretization.rectangular_mesh.size;

    scalar_value = std::min(scalar_value, number(1.0));
  }

  void
  compute_rhs(FieldContainer<dim, degree, number> &variable_list,
              const SimulationTimer               &sim_timer,
              unsigned int                         solve_block_id) const override
  {
    SystemContainer<dim, degree, number> sys_container(sys, get_user_inputs());
    if (solve_block_id == ParaboloidSystem::explicit_block_id)
      {
        sys_container.initialize_fields_explicit(variable_list);

        sys_container.calculate_sum_sq_eta();
        sys_container.calculate_h();
        sys_container.calculate_dhdeta();
        sys_container.calculate_local_mobility();
        sys_container.calculate_dmudt();

        sys_container.submit_fields_explicit(variable_list, sim_timer.get_timestep());
      }
    else if (solve_block_id == ParaboloidSystem::detadt_block_id)
      {
        sys_container.initialize_fields_nonexplicit(variable_list);

        sys_container.calculate_omega_phase();
        sys_container.calculate_sum_sq_eta();
        sys_container.calculate_h();
        sys_container.calculate_dhdeta();
        sys_container.calculate_detadt();

        sys_container.submit_fields_aux(variable_list);
      }
    else if (solve_block_id == ParaboloidSystem::pp_block_id)
      {
        sys_container.initialize_fields_postprocess(variable_list);
        sys_container.calculate_sum_sq_eta();
        sys_container.calculate_h();
        sys_container.submit_fields_postprocess(variable_list);
      }
  }

  number m_well;
  number kappa;
};

PRISMS_PF_END_NAMESPACE
