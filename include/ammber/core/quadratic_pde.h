#pragma once

#include <ammber/core/quadratic_equations.h>

#include <prismspf/core/pde_operator_base.h>

#include <deal.II/base/exceptions.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number, unsigned int num_comps>
class QuadraticGrandPotentialPDE : public PDEOperatorBase<dim, degree, number>
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
  QuadraticSystem<num_comps> sys;

  /**
   * @brief Constructor.
   */
  QuadraticGrandPotentialPDE(const UserInputParameters<dim>   &_user_inputs,
                             PhaseFieldTools<dim>             &_pf_tools,
                             const QuadraticSystem<num_comps> &_sys)
    : PDEOperatorBase<dim, degree, number>(_user_inputs, _pf_tools)
    , sys(_sys)
  {}

  void
  compute_primary_rhs(FieldContainer<dim, degree, number>                &variable_list,
                      const SimulationTimer                              &sim_timer,
                      QuadraticEquations<dim, degree, number, num_comps> &sys_container) const
  {
    sys_container.initialize_fields_explicit(variable_list);

    sys_container.calculate_c_phase();
    sys_container.calculate_sum_sq_eta();
    sys_container.calculate_h();
    sys_container.calculate_dhdeta();
    sys_container.calculate_local_mobility();
    sys_container.calculate_dmudt();

    sys_container.submit_fields_explicit(variable_list, sim_timer.get_timestep());
  }

  void
  compute_detadt_rhs(FieldContainer<dim, degree, number>                &variable_list,
                     const SimulationTimer                              &sim_timer,
                     QuadraticEquations<dim, degree, number, num_comps> &sys_container) const
  {
    sys_container.initialize_fields_aux(variable_list);

    sys_container.calculate_c_phase();
    sys_container.calculate_omega_phase();
    sys_container.calculate_sum_sq_eta();
    sys_container.calculate_h();
    sys_container.calculate_dhdeta();
    sys_container.calculate_detadt();

    sys_container.submit_fields_aux(variable_list);
  }

  void
  compute_postprocess_rhs(FieldContainer<dim, degree, number>                &variable_list,
                          const SimulationTimer                              &sim_timer,
                          QuadraticEquations<dim, degree, number, num_comps> &sys_container) const
  {
    sys_container.initialize_fields_postprocess(variable_list);

    sys_container.calculate_c_phase();
    sys_container.calculate_omega_phase();
    sys_container.calculate_sum_sq_eta();
    sys_container.calculate_h();

    sys_container.submit_fields_postprocess(variable_list);
  }

  void
  compute_rhs(FieldContainer<dim, degree, number> &variable_list,
              const SimulationTimer               &sim_timer,
              unsigned int                         solve_block_id) const override
  {
    QuadraticEquations<dim, degree, number, num_comps> sys_container(sys);
    if (solve_block_id == QuadraticSystem<num_comps>::explicit_block_id)
      {
        compute_primary_rhs(variable_list, sim_timer, sys_container);
      }
    else if (solve_block_id == QuadraticSystem<num_comps>::detadt_block_id)
      {
        compute_detadt_rhs(variable_list, sim_timer, sys_container);
      }
    else if (solve_block_id == QuadraticSystem<num_comps>::pp_block_id)
      {
        compute_postprocess_rhs(variable_list, sim_timer, sys_container);
      }
  }

  /**
   * @brief Given the order parameter values provided in `eta_0`,
   * submit the initial values and compositions to PRISMS-PF
   */
  void
  submit_ic_from_fields(const std::vector<number> &eta_0,
                        const unsigned int        &index,
                        number                    &scalar_value) const
  {
    Assert(eta_0.size() == sys.num_ops(), dealii::ExcMessage("Order parameter vector is incorrect size!"));
    const double sum_sq_eta = sum_sq(eta_0) + 1e-8;
    for (uint comp_index = 0; comp_index < sys.num_comps(); comp_index++)
      {
        unsigned int var_index = sys.mu_base() + comp_index;
        if (index == var_index)
          {
            double mu0          = 0.0;
            double k_inv_interp = 0.0;
            for (unsigned int op_index = 0; op_index < sys.num_ops(); op_index++)
              {
                unsigned int phase_index     = sys.order_params[op_index];
                auto        &phase_comp_info = sys.phases.at(phase_index).comps.at(comp_index);
                mu0 += (eta_0[op_index] * eta_0[op_index] / sum_sq_eta) *
                       (phase_comp_info.c0 - phase_comp_info.c_min);
                k_inv_interp += (eta_0[op_index] * eta_0[op_index] / sum_sq_eta) / phase_comp_info.k_well;
              }
            mu0 /= k_inv_interp;
            scalar_value = mu0;
            return;
          }
      }
    for (unsigned int op_index = 0; op_index < sys.num_ops(); op_index++)
      {
        if (index == sys.eta_base() + op_index)
          {
            scalar_value = eta_0[op_index];
            return;
          }
      }
  }

  template <typename vectorType>
  auto
  sum_sq(const vectorType &vec) const
  {
    decltype(vec[0] * vec[0]) sum = 0.0;
    for (unsigned int i = 0; i < vec.size(); i++)

      {
        const auto &val = vec[i];
        sum += val * val;
      }
    return sum;
  }
};

PRISMS_PF_END_NAMESPACE
