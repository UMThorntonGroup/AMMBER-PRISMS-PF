#ifndef SYSTEMCONTAINER_H
#define SYSTEMCONTAINER_H

#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

#include "ParaboloidSystem.h"

#include <FieldContainer.h>
#include <core/userInputParameters.h>
#include <core/variableContainer.h>
#include <map>
#include <memory>
#include <string>

/**
 * @brief Class for solving the PDE (equations.cc and postprocess.cc)
 * @tparam dim The dimension of the problem
 * @tparam degree The degree of the finite element
 */
template <int dim, int degree>
class SystemContainer
{
public:
  using scalarValue     = dealii::VectorizedArray<double>;
  using scalarGrad      = dealii::Tensor<1, dim, scalarValue>;
  using scalarField     = Dual<scalarValue, dim>;
  using scalarVariation = Variation<scalarValue, dim>;

  /**
   * @brief Data structure to hold the phase data
   */
  struct PhaseData
  {
    scalarField omega;
    scalarField h;
  };

  /**
   * @brief Data structure to hold the composition data
   */
  struct CompData
  {
    scalarField     mu;
    scalarVariation dmudt;
    scalarValue     M;
  };

  /**
   * @brief Data structure to hold the order parameter data
   */
  struct OPData
  {
    scalarField              eta;
    scalarVariation          detadt;
    scalarValue              detadt_field = dealii::make_vectorized_array(0.);
    std::vector<scalarField> dhdeta;
  };

  /**
   * @brief Pointer to the system parameters
   */
  std::shared_ptr<const ParaboloidSystem> isoSys;
  /**
   * @brief Pointer to the PRISMS-PF parameters
   */
  std::shared_ptr<const userInputParameters<dim>> userInputs;

  /**
   * @brief Values associated with each phase
   */
  std::vector<PhaseData> phase_data;
  /**
   * @brief Values associated with each component
   */
  std::vector<CompData> comp_data;
  /**
   * @brief Values associated with each order parameter
   */
  std::vector<std::pair<uint, OPData>> op_data;
  /**
   * @brief Sum of squares of the order parameters
   */
  scalarField sum_sq_eta;

  /**
   * Constructor
   */
  SystemContainer(const ParaboloidSystem &sys, const userInputParameters<dim> &inputs)
    : isoSys(std::make_shared<const ParaboloidSystem>(sys))
    , userInputs(std::make_shared<const userInputParameters<dim>>(inputs))
    , phase_data(std::vector<PhaseData>(isoSys->phases.size()))
    , comp_data(std::vector<CompData>(isoSys->comp_names.size()))
    , op_data({})
    , sum_sq_eta({})
  {}

  ~SystemContainer()
  {}

  /**
   * @brief Initialize the fields for the PDE
   * @param variable_list The variable list
   * @param var_index The starting index for the block of fields
   */
  void
  initialize_fields_explicit(
    const variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
    uint                                                                  &var_index)
  {
    op_data.clear();
    op_data.reserve(isoSys->order_params.size());
    for (uint comp_index = 0; comp_index < isoSys->comp_names.size(); comp_index++)
      {
        comp_data[comp_index].mu.val  = variable_list.get_scalar_value(var_index);
        comp_data[comp_index].mu.grad = variable_list.get_scalar_gradient(var_index);
        var_index++;
      }
    for (const auto &phase_index : isoSys->order_params)
      {
        OPData op;
        op.eta.val  = variable_list.get_scalar_value(var_index);
        op.eta.grad = variable_list.get_scalar_gradient(var_index);
        op.detadt_field =
          variable_list.get_scalar_value(var_index + isoSys->order_params.size());
        op.dhdeta.resize(isoSys->phases.size());
        op_data.push_back({phase_index, op});
        var_index++;
      }
  }

  /**
   * @brief Initialize the fields needed for the postprocess
   * @param variable_list The variable list
   * @param var_index The starting index for the block of fields
   */
  void
  initialize_fields_postprocess(
    const variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
    uint                                                                  &var_index)
  {
    op_data.clear();
    op_data.reserve(isoSys->order_params.size());
    for (uint comp_index = 0; comp_index < isoSys->comp_names.size(); comp_index++)
      {
        comp_data[comp_index].mu.val = variable_list.get_scalar_value(var_index);
        var_index++;
      }
    for (const auto &phase_index : isoSys->order_params)
      {
        OPData op;
        op.eta.val = variable_list.get_scalar_value(var_index);
        op.dhdeta.resize(isoSys->phases.size());
        op_data.push_back({phase_index, op});
        var_index++;
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
        const ParaboloidSystem::Phase &phase_info = isoSys->phases[phase_index];
        PhaseData                     &phase      = phase_data[phase_index];
        phase.omega.val                           = phase_info.f_min;
        for (uint comp_index = 0; comp_index < comp_data.size(); comp_index++)
          {
            const CompData                        &comp = comp_data.at(comp_index);
            const ParaboloidSystem::PhaseCompInfo &comp_info =
              phase_info.comps.at(comp_index);
            phase.omega +=
              -comp.mu * comp.mu / (2.0 * comp_info.k_well) - comp.mu * comp_info.c_min;
          }
      }
  }

  /**
   * @brief Calculate the sum of squares of the order parameters
   */
  void
  calculate_sum_sq_eta()
  {
    sum_sq_eta.val = dealii::make_vectorized_array(0.);
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
            scalarField &dhdeta = op.dhdeta[beta_index];
            dhdeta.val          = dealii::make_vectorized_array(0.);
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
    scalarValue L               = scalarValue(0.0);
    scalarValue sum_pair_sq_eta = scalarValue(0.0);
    // scalarValue sq_sum_sq_eta   = sum_sq_eta.val * sum_sq_eta.val;
    for (const auto &[alpha_index, op1] : op_data)
      {
        for (const auto &[beta_index, op2] : op_data)
          {
            if (&op1 == &op2)
              {
                continue; // Skip self-interaction
              }
            double mu_ab =
              2.0 *
              (isoSys->phases[alpha_index].mu_int * isoSys->phases[beta_index].mu_int) /
              (isoSys->phases[alpha_index].mu_int + isoSys->phases[beta_index].mu_int);
            double L_ab = 4.0 * mu_ab / isoSys->l_int / 3.0;
            L += L_ab * (op1.eta.val * op1.eta.val + op2.eta.val * op2.eta.val);
            sum_pair_sq_eta += op1.eta.val * op1.eta.val + op2.eta.val * op2.eta.val;
          }
      }
    L /= 2.0 * sum_pair_sq_eta + 1.0e-8;
    for (auto &[alpha_index, op] : op_data)
      {
        const ParaboloidSystem::Phase &phase_info = isoSys->phases.at(alpha_index);

        double m     = 6.00 * phase_info.sigma / isoSys->l_int;
        double kappa = 0.75 * phase_info.sigma * isoSys->l_int;

        // Interface term
        scalarVariation interface_term;
        interface_term.val =
          m * (op.eta.val * op.eta.val * op.eta.val - op.eta.val +
               2. * 1.5 * op.eta.val * (sum_sq_eta.val - op.eta.val * op.eta.val));
        interface_term.vec = -kappa * op.eta.grad;

        // Chemical term
        // This is a variation, but has no vector term.
        scalarValue chemical_term = dealii::make_vectorized_array(0.);
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
    for (uint comp_index = 0; comp_index < comp_data.size(); comp_index++)
      {
        CompData &comp = comp_data[comp_index];
        comp.M         = dealii::make_vectorized_array(0.);
        for (uint phase_index = 0; phase_index < phase_data.size(); phase_index++)
          {
            PhaseData &phase = phase_data[phase_index];
            comp.M += isoSys->phases.at(phase_index).D * phase.h.val /
                      (isoSys->phases.at(phase_index).comps.at(comp_index).k_well);
          }
      }
  }

  /**
   * @brief Calculate the time evolution of the chemical potential
   * @details Calculates the evolution of the composition, then converts to chemical
   * potential through the susceptibility, chi_AA
   */
  void
  calculate_dmudt()
  {
    for (uint comp_index = 0; comp_index < comp_data.size(); comp_index++)
      {
        CompData &comp = comp_data[comp_index];

        // Calculate the susceptibility
        scalarField chi_AA;
        for (uint phase_index = 0; phase_index < phase_data.size(); phase_index++)
          {
            PhaseData                             &phase = phase_data[phase_index];
            const ParaboloidSystem::PhaseCompInfo &comp_info =
              isoSys->phases.at(phase_index).comps.at(comp_index);
            chi_AA += phase.h / comp_info.k_well;
          }

        // Flux term
        comp.dmudt.val = dealii::make_vectorized_array(0.);
        comp.dmudt.vec = -comp.M * -comp.mu.grad;

        // Partitioning term
        for (auto &[phase_index, op] : op_data)
          {
            scalarField dcdeta_sum;
            for (uint beta_index = 0; beta_index < phase_data.size(); beta_index++)
              {
                auto &comp_info = isoSys->phases.at(beta_index).comps.at(comp_index);
                dcdeta_sum += op.dhdeta.at(beta_index) *
                              (comp.mu / comp_info.k_well + comp_info.c_min);
              }
            comp.dmudt -= dcdeta_sum.val * op.detadt_field;
          }

        // Convert from dcdt to dmudt
        comp.dmudt = (1.0 / chi_AA) * comp.dmudt;
      }
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
   * @param var_index The starting index for the block of fields
   */
  void
  submit_fields(
    variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
    uint                                                            &var_index)
  {
    for (uint comp_index = 0; comp_index < comp_data.size(); comp_index++)
      {
        CompData &comp = comp_data[comp_index];
        variable_list.set_scalar_value_term_RHS(var_index,
                                                comp.mu.val +
                                                  comp.dmudt.val * userInputs->dtValue);
        variable_list.set_scalar_gradient_term_RHS(var_index,
                                                   -comp.dmudt.vec * userInputs->dtValue);
        var_index++;
      }
    for (auto &[phase_index, op] : op_data)
      {
        variable_list.set_scalar_value_term_RHS(var_index,
                                                op.eta.val +
                                                  op.detadt_field * userInputs->dtValue);
        variable_list.set_scalar_gradient_term_RHS(var_index, -op.detadt.vec * 0.0);
        var_index++;
      }
  }

  /**
   * @brief Submit the fields to PRISMS-PF
   * @param variable_list The PRISMS-PF variable list
   * @param var_index The starting index for the block of fields
   */
  void
  submit_aux_fields(
    variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
    uint                                                            &var_index)
  {
    for (uint comp_index = 0; comp_index < comp_data.size(); comp_index++)
      {
        var_index++;
      }
    for (auto &[phase_index, op] : op_data)
      {
        var_index++;
      }
    for (auto &[phase_index, op] : op_data)
      {
        variable_list.set_scalar_value_term_RHS(var_index, op.detadt.val);
        variable_list.set_scalar_gradient_term_RHS(var_index, -op.detadt.vec);
        var_index++;
      }
  }

  /**
   * @brief Submit the post-processed fields to PRISMS-PF
   * @param pp_variable_list The PRISMS-PF variable list
   * @param pp_index The starting index for the block of fields
   */
  void
  submit_pp_fields(
    variableContainer<dim, degree, dealii::VectorizedArray<double>> &pp_variable_list,
    uint                                                            &pp_index)
  {
    for (uint comp_index = 0; comp_index < comp_data.size(); comp_index++)
      {
        CompData   &comp = comp_data[comp_index];
        scalarValue c    = dealii::make_vectorized_array(0.);
        for (uint phase_index = 0; phase_index < phase_data.size(); phase_index++)
          {
            const PhaseData                       &phase = phase_data.at(phase_index);
            const ParaboloidSystem::PhaseCompInfo &comp_info =
              isoSys->phases.at(phase_index).comps.at(comp_index);
            c += phase.h.val * comp_info.c_min;
            c += phase.h.val * comp.mu.val / comp_info.k_well;
          }
        pp_variable_list.set_scalar_value_term_RHS(pp_index, c);
        pp_index++;
      }
  }
};

#endif