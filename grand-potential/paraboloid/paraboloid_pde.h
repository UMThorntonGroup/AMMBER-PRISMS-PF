#include "system_equations.h"

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

  /**
   * @brief Object containing the thermodynamic and kinetic parameters
   */
  ParaboloidSystem sys;

  /**
   * @brief Constructor.
   */
  CustomPDE(const UserInputParameters<dim> &_user_inputs,
            PhaseFieldTools<dim>           &_pf_tools,
            const ParaboloidSystem         &_sys)
    : PDEOperatorBase<dim, degree, number>(_user_inputs, _pf_tools)
    , sys(_sys)
  {}

private:
  void
  set_initial_condition([[maybe_unused]] const unsigned int       &index,
                        [[maybe_unused]] const unsigned int       &component,
                        [[maybe_unused]] const dealii::Point<dim> &point,
                        [[maybe_unused]] number                   &scalar_value,
                        [[maybe_unused]] number                   &vector_component_value) const override
  {
    const dealii::Tensor<1, dim> &mesh_size = get_user_inputs().spatial_discretization.rectangular_mesh.size;

    scalar_value = 0.0;
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
        sys_container.initialize_fields_aux(variable_list);

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

  /**
   * @brief Prints the grand potential densities (nondimensionalized) of each phase at
   * its initial conditions
   */
  void
  print_initial_energies()
  {
    /**
     * @brief Map of the initial grand potential densities for each phase (used for
     * printing)
     */
    /* std::map<std::string, double> initial_omega;
    std::cout << "Initial omega free energies:\n";
    for (uint phase_index = 0; phase_index < sys.phases.size(); phase_index++)
      {
        const ParaboloidSystem::Phase       &phase = sys.phases.at(phase_index);
        SystemContainer<dim, degree, number> sys_for_print(sys);
        sys_for_print.op_data.push_back({
          phase_index,
          {
            {ScalarValue(1.0), {}}, // eta
            {ScalarValue(0.0), {}}, // detadt
            ScalarValue(0.0),       // detadt_field
            {}                      // dhdeta
          }
        });

        for (uint comp_index = 0; comp_index < sys.comp_names.size(); comp_index++)
          {
            const ParaboloidSystem::PhaseCompInfo &comp_info = phase.comps.at(comp_index);
            double                                 mu0 = comp_info.k_well * (comp_info.x0 - comp_info.c_min);
            sys_for_print.comp_data[comp_index].mu.val = ScalarValue(mu0);
          }
        sys_for_print.calculate_omega_phase();

        std::cout << phase.name << ":\n"
                  << "Omega:\t" << sys_for_print.phase_data[phase_index].omega.val[0] << "\n";
        initial_omega[phase.name] = sys_for_print.phase_data[phase_index].omega.val[0];
        for (uint comp_index = 0; comp_index < sys.comp_names.size(); comp_index++)
          {
            const ParaboloidSystem::PhaseCompInfo &comp = phase.comps.at(comp_index);
            std::cout << "mu_" << comp.name << ":\t" << sys_for_print.comp_data[comp_index].mu.val[0] << "\n";
          }
        std::cout << "\n";
      } */
  }

  /**
   * @brief Function to print the properties of the interface between two phases at
   * initial conditions
   */
  /*   void
    print_interface_properties()
    {
      for (const auto &alpha : sys.phases)
        {
          for (const auto &beta : sys.phases)
            {
              std::cout << "Properties of the interface between " << alpha.name << " and " << beta.name
                        << ":\n";
              double delta_g = initial_omega[alpha.name] - initial_omega[beta.name];
              double sigma   = 0.5 * (alpha.sigma + beta.sigma);
              double D       = 0.5 * (alpha.D * beta.D) / (alpha.D + beta.D);
              double mu_int  = 0.5 * (alpha.mu_int * beta.mu_int) / (alpha.mu_int + beta.mu_int);
              std::cout << "The value of the dimensionless number delta_g/(sigma/l_int) is "
                        << delta_g * sys.l_int / sigma << "\n";
              std::cout << "The value of the dimensionless number delta_g*mu_int*l_int/D is "
                        << delta_g * mu_int * sys.l_int / D << "\n";
              std::cout << "The value of the dimensionless number sigma*mu_int/D is " << sigma * mu_int / D
                        << "\n\n";
            }
        }
    } */
};

PRISMS_PF_END_NAMESPACE
