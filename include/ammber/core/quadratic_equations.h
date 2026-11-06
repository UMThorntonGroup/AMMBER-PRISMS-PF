#pragma once

#include <ammber/core/quadratic_system.h>
#include <ammber/utilities/auto_diff.h>
#include <ammber/utilities/mat_inversion.h>

#include <prismspf/core/field_container.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

using namespace prismspf;

/**
 * @brief Class for solving the PDE (equations.cc and postprocess.cc)
 * @tparam dim The dimension of the problem
 * @tparam degree The degree of the finite element
 * @tparam num_comps The number of components
 */
template <int dim, int degree, typename number, unsigned num_comps>
class QuadraticEquations
{
public:
  using ScalarValue     = dealii::VectorizedArray<number>;
  using ScalarGrad      = dealii::Tensor<1, dim, ScalarValue>;
  using ScalarField     = Dual<ScalarValue, dim>;
  using ScalarVariation = Variation<ScalarValue, dim>;
  using VarList         = FieldContainer<dim, degree, number>;
  template <typename T>
  using CompVector = dealii::Tensor<1, num_comps, T>;
  template <typename T>
  using CompMatrix = dealii::Tensor<2, num_comps, T>;
  template <typename T>
  using CompSymmet = dealii::SymmetricTensor<2, num_comps, T>;

  /**
   * @brief Data structure to hold the phase data
   */
  struct PhaseData
  {
    ScalarField             omega;
    ScalarField             h;
    CompVector<ScalarField> c_phase;
  };

  /**
   * @brief Data structure to hold the order parameter data
   */
  struct OPData
  {
    ScalarField              eta;
    ScalarVariation          detadt;
    ScalarValue              detadt_field = 0.0;
    std::vector<ScalarField> dhdeta;
  };

  /**
   * @brief Pointer to the system parameters
   */
  const QuadraticSystem<num_comps> *sys_ptr;

  /**
   * @brief Values associated with each phase
   */
  std::vector<PhaseData> phase_data;
  /**
   * @brief Values associated with each component
   */
  CompVector<ScalarField>     mu;
  CompVector<ScalarVariation> dmudt;
  CompMatrix<ScalarValue>     M;

  /**
   * @brief Values associated with each order parameter
   */
  std::vector<std::pair<uint, OPData>> op_data;

  /**
   * @brief Sum of squares of the order parameters
   */
  ScalarField sum_sq_eta;

  /**
   * Constructor
   */
  explicit QuadraticEquations(const QuadraticSystem<num_comps> &_sys)
    : sys_ptr(&_sys)
    , phase_data(sys().phases.size())
    , op_data(sys().num_ops())
    , sum_sq_eta({})
  {}

  ~QuadraticEquations()
  {}

  const QuadraticSystem<num_comps> &
  sys() const
  {
    return *sys_ptr;
  }

  /**
   * @brief Initialize the fields for the PDE
   * @param variable_list The variable list
   */
  void
  initialize_fields_explicit(const VarList &variable_list)
  {
    for (uint comp_index = 0; comp_index < num_comps; comp_index++)
      {
        mu[comp_index].val = variable_list.template get_value<Scalar, OldOne>(sys().mu_base() + comp_index);
        mu[comp_index].grad =
          variable_list.template get_gradient<Scalar, OldOne>(sys().mu_base() + comp_index);
      }
    for (unsigned int op_index = 0; op_index < sys().num_ops(); op_index++)
      {
        const uint phase_index  = sys().order_params[op_index];
        OPData    &op           = op_data[op_index].second;
        op_data[op_index].first = phase_index;
        op.eta.val      = variable_list.template get_value<Scalar, OldOne>(sys().eta_base() + op_index);
        op.eta.grad     = variable_list.template get_gradient<Scalar, OldOne>(sys().eta_base() + op_index);
        op.detadt_field = variable_list.template get_value<Scalar, OldOne>(sys().detadt_base() + op_index);
        op.dhdeta.resize(sys().phases.size());
      }
  }

  /**
   * @brief Initialize the fields for the PDE
   * @param variable_list The variable list
   */
  void
  initialize_fields_aux(const VarList &variable_list)
  {
    for (uint comp_index = 0; comp_index < num_comps; comp_index++)
      {
        mu[comp_index].val = variable_list.template get_value<Scalar, Current>(sys().mu_base() + comp_index);
        mu[comp_index].grad =
          variable_list.template get_gradient<Scalar, Current>(sys().mu_base() + comp_index);
      }
    for (unsigned int op_index = 0; op_index < sys().num_ops(); op_index++)
      {
        const uint phase_index  = sys().order_params[op_index];
        OPData    &op           = op_data[op_index].second;
        op_data[op_index].first = phase_index;
        op.eta.val  = variable_list.template get_value<Scalar, Current>(sys().eta_base() + op_index);
        op.eta.grad = variable_list.template get_gradient<Scalar, Current>(sys().eta_base() + op_index);
        op.dhdeta.resize(sys().phases.size());
      }
  }

  /**
   * @brief Initialize the fields needed for the postprocess
   * @param variable_list The variable list
   * @param var_index The starting index for the block of fields
   */
  void
  initialize_fields_postprocess(const VarList &variable_list)
  {
    for (uint comp_index = 0; comp_index < num_comps; comp_index++)
      {
        mu[comp_index].val = variable_list.template get_value<Scalar, Current>(sys().mu_base() + comp_index);
      }
    for (unsigned int op_index = 0; op_index < sys().num_ops(); op_index++)
      {
        const uint phase_index  = sys().order_params[op_index];
        op_data[op_index].first = phase_index;
        OPData &op              = op_data[op_index].second;
        op.eta.val = variable_list.template get_value<Scalar, Current>(sys().eta_base() + op_index);
        op.dhdeta.resize(sys().phases.size());
      }
  }

  /**
   * @brief Calculate the grand potential density for each phase
   */
  void
  calculate_c_phase()
  {
    for (uint phase_index = 0; phase_index < phase_data.size(); phase_index++)
      {
        const typename QuadraticSystem<num_comps>::Phase &phase_info = sys().phases[phase_index];
        PhaseData                                        &phase      = phase_data[phase_index];
        phase.omega.val                                              = phase_info.C;
        auto delta_mu                                                = mu - phase_info.B;
        phase.c_phase = phase_info.c_ref + phase_info.suscept * delta_mu;
      }
  }

  /**
   * @brief Calculate the grand potential density for each phase
   */
  void
  calculate_omega_phase()
  {
    for (uint phase_index = 0; phase_index < phase_data.size(); phase_index++)
      {
        const typename QuadraticSystem<num_comps>::Phase &phase_info = sys().phases[phase_index];
        PhaseData                                        &phase      = phase_data[phase_index];
        phase.omega.val                                              = phase_info.C;
        auto delta_mu                                                = mu - phase_info.B;
        auto delta_c                                                 = phase_info.suscept * delta_mu;
        phase.omega += -0.5 * delta_mu * delta_c - mu * phase.c_phase;
      }
  }

  /**
   * @brief Calculate the sum of squares of the order parameters
   */
  void
  calculate_sum_sq_eta()
  {
    sum_sq_eta.val = 0.0;
    for (const auto &[phase_index, op] : op_data)
      {
        sum_sq_eta += op.eta * op.eta;
      }
  }

  /**
   * @brief Calculate the phase fraction, h, for each phase
   * @details h_i = eta_i^2 / sum(eta^2)
   */
  void
  calculate_h()
  {
    for (auto &[phase_index, op] : op_data)
      {
        phase_data[phase_index].h += op.eta * op.eta;
      }
    for (auto &phase : phase_data)
      {
        phase.h /= sum_sq_eta;
      }
  }

  /**
   * @brief Calculate the derivative of the phase fraction, h, with respect to eta
   */
  void
  calculate_dhdeta()
  {
    for (auto &[alpha_index, op] : op_data)
      {
        for (uint beta_index = 0; beta_index < op.dhdeta.size(); beta_index++)
          {
            PhaseData   &beta   = phase_data[beta_index];
            ScalarField &dhdeta = op.dhdeta[beta_index];
            dhdeta.val          = 0.0;
            if (alpha_index == beta_index)
              {
                dhdeta += 2.0 * op.eta;
              }
            dhdeta -= 2.0 * op.eta * beta.h;
            dhdeta /= sum_sq_eta;
          }
      }
  }

  /**
   * @brief Calculate the time evolution of the order parameters
   */
  void
  calculate_detadt()
  {
    // Calculate the local interface mobility
    ScalarValue m(0.0);
    ScalarValue kappa(0.0);
    ScalarValue L(0.0);
    ScalarValue sum_pair_sq_eta = ScalarValue(0.0);
    // ScalarValue sq_sum_sq_eta   = sum_sq_eta.val * sum_sq_eta.val;
    for (const auto &[alpha_index, op1] : op_data)
      {
        for (const auto &[beta_index, op2] : op_data)
          {
            if (&op1 == &op2)
              {
                continue; // Skip self-interaction
              }
            double mu_ab = 2.0 * (sys().phases[alpha_index].mu_int * sys().phases[beta_index].mu_int) /
                           (sys().phases[alpha_index].mu_int + sys().phases[beta_index].mu_int);
            double L_ab = 4.0 * mu_ab / sys().l_int / 3.0;
            double sigma_ab =
              0.5 * (sys().sigma(alpha_index, beta_index) + sys().sigma(beta_index, alpha_index));
            double      m_ab        = 6.00 * sigma_ab / sys().l_int;
            double      kappa_ab    = 0.75 * sigma_ab * sys().l_int;
            ScalarValue pair_sq_eta = op1.eta.val * op1.eta.val * op2.eta.val * op2.eta.val;
            L += L_ab * pair_sq_eta;
            m += m_ab * pair_sq_eta;
            kappa += kappa_ab * pair_sq_eta;
            sum_pair_sq_eta += pair_sq_eta;
          }
      }
    L /= sum_pair_sq_eta + 1.0e-12;
    m /= sum_pair_sq_eta + 1.0e-12;
    kappa /= sum_pair_sq_eta + 1.0e-12;
    for (auto &[alpha_index, op] : op_data)
      {
        const typename QuadraticSystem<num_comps>::Phase &phase_info = sys().phases.at(alpha_index);

        // Interface term
        ScalarVariation interface_term;
        interface_term.val = m * (op.eta.val * op.eta.val * op.eta.val - op.eta.val +
                                  2. * 1.5 * op.eta.val * (sum_sq_eta.val - op.eta.val * op.eta.val));
        interface_term.vec = -kappa * op.eta.grad;

        // Chemical term
        // This is a variation, but has no vector term.
        ScalarValue chemical_term(0.0);
        for (uint beta_index = 0; beta_index < phase_data.size(); beta_index++)
          {
            const PhaseData &beta = phase_data.at(beta_index);
            chemical_term += beta.omega.val * op.dhdeta.at(beta_index).val;
          }
        op.detadt = -L * (interface_term + chemical_term);
      }
  }

  /**
   * @brief Calculate the local chemical mobility
   */
  void
  calculate_local_mobility()
  {
    M.clear();
    for (uint phase_index = 0; phase_index < phase_data.size(); phase_index++)
      {
        PhaseData                                        &phase      = phase_data[phase_index];
        const typename QuadraticSystem<num_comps>::Phase &phase_info = sys().phases[phase_index];
        const CompMatrix<double>                         &D          = phase_info.D;
        const CompMatrix<double>                         &suscept    = phase_info.suscept;
        M += phase.h.val * phase_info.D * phase_info.suscept;
      }
  }

  /**
   * @brief Calculate the time evolution of the chemical potential
   * @details Calculates the evolution of the composition, then converts to chemical
   * potential through the susceptibility, local_suscept
   */
  void
  calculate_dmudt()
  {
    CompSymmet<ScalarField> local_suscept;
    for (uint phase_index = 0; phase_index < phase_data.size(); phase_index++)
      {
        PhaseData                                        &phase      = phase_data[phase_index];
        const typename QuadraticSystem<num_comps>::Phase &phase_info = sys().phases[phase_index];
        local_suscept += phase.h * dealii::symmetrize(phase_info.suscept);
      }
    // Flux term
    for (uint comp_i = 0; comp_i < num_comps; comp_i++)
      {
        dmudt[comp_i].val = 0.0;
        for (uint comp_j = 0; comp_j < num_comps; comp_j++)
          {
            dmudt[comp_i].vec += M * mu[comp_j].grad;
          }
      }
    // Partitioning term
    for (auto &[phase_index, op] : op_data)
      {
        CompVector<ScalarField> dcdeta;
        for (uint beta_index = 0; beta_index < phase_data.size(); beta_index++)
          {
            dcdeta += op.dhdeta.at(beta_index) * (phase_data[beta_index].c_phase);
          }
        for (uint comp_index = 0; comp_index < num_comps; comp_index++)
          {
            dmudt[comp_index] -= op.detadt_field * dcdeta[comp_index].val;
          }
      }
    // Convert from dcdt to dmudt
    CompMatrix<ScalarField> suscept_inv = inverse_from_cholesky(local_suscept);
    dmudt                               = suscept_inv * dmudt;
  }

  /**
   * @brief Calculate the information needed to solve the evolution equations in the
   * proper order
   */
  void
  calculate_locals()
  {
    calculate_omega_phase();
    calculate_sum_sq_eta();
    calculate_h();
    calculate_dhdeta();
    calculate_local_mobility();
  }

  /**
   * @brief Submit the fields to PRISMS-PF
   * @param variable_list The PRISMS-PF variable list
   * @param dt The timestep
   */
  void
  submit_fields_explicit(VarList &variable_list, double dt)
  {
    for (uint comp_index = 0; comp_index < num_comps; comp_index++)
      {
        variable_list.set_value_term(sys().mu_base() + comp_index,
                                     mu[comp_index].val + dmudt[comp_index].val * dt);
        variable_list.set_gradient_term(sys().mu_base() + comp_index, -mu[comp_index].grad * dt);
      }
    for (uint op_index = 0; op_index < sys().num_ops(); op_index++)
      {
        const auto &[phase_index, op] = op_data[op_index];
        variable_list.set_value_term(sys().eta_base() + op_index, op.eta.val + op.detadt_field * dt);
      }
  }

  /**
   * @brief Submit the fields to PRISMS-PF
   * @param variable_list The PRISMS-PF variable list
   * @param var_index The starting index for the block of fields
   */
  void
  submit_fields_aux(VarList &variable_list)
  {
    for (uint op_index = 0; op_index < sys().num_ops(); op_index++)
      {
        const auto &[phase_index, op] = op_data[op_index];
        variable_list.set_value_term(sys().detadt_base() + op_index, op.detadt.val);
        variable_list.set_gradient_term(sys().detadt_base() + op_index, -op.detadt.vec);
      }
  }

  /**
   * @brief Submit the post-processed fields to PRISMS-PF
   * @param variable_list The PRISMS-PF variable list
   */
  void
  submit_fields_postprocess(VarList &variable_list)
  {
    for (uint comp_index = 0; comp_index < num_comps; comp_index++)
      {
        ScalarValue c(0.0);
        for (uint phase_index = 0; phase_index < sys().num_phases(); phase_index++)
          {
            const PhaseData &phase = phase_data.at(phase_index);
            c += phase.h.val * phase.c_phase[comp_index].val;
          }
        variable_list.set_value_term(sys().c_tot_base() + comp_index, c);
      }
    for (uint phase_index = 0; phase_index < sys().num_phases(); phase_index++)
      {
        const PhaseData &phase = phase_data.at(phase_index);
        variable_list.set_value_term(sys().h_base() + phase_index, phase.h.val);
        variable_list.set_value_term(sys().omega_base() + phase_index, phase.omega.val);
      }
  }
};
