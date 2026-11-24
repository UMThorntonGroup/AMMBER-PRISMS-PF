#pragma once

#include <ammber/core/quadratic_equations.h>

#include <prismspf/core/pde_operator_base.h>

#include <deal.II/base/exceptions.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number, typename system_type>
class QuadraticGrandPotentialPDE : public PDEOperatorBase<dim, degree, number>
{
public:
  static constexpr unsigned int num_comps = system_type::num_comps;
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
  system_type sys;

  /**
   * @brief Constructor.
   */
  QuadraticGrandPotentialPDE(const UserInputParameters<dim>   &_user_inputs,
                             PhaseFieldTools<dim>             &_pf_tools,
                             const system_type &_sys)
    : PDEOperatorBase<dim, degree, number>(_user_inputs, _pf_tools)
    , sys(_sys)
  {}

  virtual ScalarValue
  get_local_temperature(FieldContainer<dim, degree, number>                &variable_list,
                        const SimulationTimer                              &sim_timer) const
  {
    return 0.0;
  }

  void
  compute_primary_rhs(FieldContainer<dim, degree, number>                  &variable_list,
                      const SimulationTimer                                &sim_timer,
                      QuadraticEquations<dim, degree, number, system_type> &sys_container) const
  {
    ScalarValue local_temperature = get_local_temperature(variable_list, sim_timer);
    sys_container.initialize_fields_explicit(variable_list);

    sys_container.calculate_c_phase(local_temperature);
    sys_container.calculate_sum_sq_eta();
    sys_container.calculate_h();
    sys_container.calculate_dhdeta();
    sys_container.calculate_local_mobility(local_temperature);
    sys_container.calculate_dmudt(local_temperature);

    sys_container.submit_fields_explicit(variable_list, sim_timer.get_timestep());
  }

  void
  compute_detadt_rhs(FieldContainer<dim, degree, number>                  &variable_list,
                     const SimulationTimer                                &sim_timer,
                     QuadraticEquations<dim, degree, number, system_type> &sys_container) const
  {
    ScalarValue local_temperature = get_local_temperature(variable_list, sim_timer);
    sys_container.initialize_fields_aux(variable_list);

    sys_container.calculate_c_phase(local_temperature);
    sys_container.calculate_omega_phase(local_temperature);
    sys_container.calculate_sum_sq_eta();
    sys_container.calculate_h();
    sys_container.calculate_dhdeta();
    sys_container.calculate_detadt(local_temperature);

    sys_container.submit_fields_aux(variable_list);
  }

  void
  compute_postprocess_rhs(FieldContainer<dim, degree, number>                  &variable_list,
                          const SimulationTimer                                &sim_timer,
                          QuadraticEquations<dim, degree, number, system_type> &sys_container) const
  {
    ScalarValue local_temperature = get_local_temperature(variable_list, sim_timer);
    sys_container.initialize_fields_postprocess(variable_list);

    sys_container.calculate_c_phase(local_temperature);
    sys_container.calculate_omega_phase(local_temperature);
    sys_container.calculate_sum_sq_eta();
    sys_container.calculate_h();

    sys_container.submit_fields_postprocess(variable_list);
  }

  void
  compute_rhs(FieldContainer<dim, degree, number> &variable_list,
              const SimulationTimer               &sim_timer,
              unsigned int                         solve_block_id) const override
  {
    QuadraticEquations<dim, degree, number, system_type> sys_container(sys);
    if (solve_block_id == system_type::explicit_block_id)
      {
        compute_primary_rhs(variable_list, sim_timer, sys_container);
      }
    else if (solve_block_id == system_type::detadt_block_id)
      {
        compute_detadt_rhs(variable_list, sim_timer, sys_container);
      }
    else if (solve_block_id == system_type::pp_block_id)
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
                        number                    &scalar_value,
                        number                     temperature) const
  {
    Assert(eta_0.size() == sys.num_ops(), dealii::ExcMessage("Order parameter vector is incorrect size!"));
    const double sum_sq_eta = sum_sq(eta_0) + 1e-8;
    // TODO: fix
    typename system_type::template CompVector<double> mu0;
    typename system_type::template CompVector<double> B_interp;
    typename system_type::template CompVector<double> c_ref_interp;
    typename system_type::template CompVector<double> c0_interp;
    typename system_type::template CompMatrix<double> sus_interp;
    for (unsigned int op_index = 0; op_index < sys.num_ops(); op_index++)
      {
        unsigned int phase_index = sys.order_params[op_index];
        const auto  &phase_info  = sys.phases[phase_index];
        number       h           = eta_0[op_index] * eta_0[op_index] / sum_sq_eta;
        B_interp += h * phase_info.get_B(temperature);
        c_ref_interp += h * phase_info.get_c_ref(temperature);
        c0_interp += h * phase_info.c0;
        sus_interp += h * phase_info.get_suscept(temperature);
      }
    mu0 = inverse_from_cholesky(dealii::symmetrize(sus_interp)) * ((c0_interp - c_ref_interp)) + B_interp;
    for (uint comp_index = 0; comp_index < num_comps; comp_index++)
      {
        unsigned int var_index = sys.mu_base() + comp_index;
        if (index == var_index)
          {
            scalar_value = mu0[comp_index];
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
