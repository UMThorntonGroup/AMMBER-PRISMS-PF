// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"
#include "system_equations.h"

#include <prismspf/config.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>
#include <prismspf/core/variable_container.h>

PRISMS_PF_BEGIN_NAMESPACE

void
CustomAttributeLoader::load_variable_attributes()
{
  // Get names for mu fields
  std::vector<std::string> mu_names;
  std::vector<std::string> grad_mu_names;
  prisms::ConditionalOStreams::pout_base() << "Comp names: ";
  for (const auto &comp_name : system.comp_names)
    {
      std::string var_name = "mu_" + comp_name;
      mu_names.push_back(var_name);
      grad_mu_names.push_back("grad(" + var_name + ")");
      prisms::ConditionalOStreams::pout_base() << var_name << " ";
    }
  prisms::ConditionalOStreams::pout_base() << "\n"
                                           << "Order parameter names: ";
  // Get names for order parameter fields
  std::vector<std::string> op_names = system.get_order_parameter_names();
  std::vector<std::string> grad_op_names;
  for (const auto &op_name : op_names)
    {
      grad_op_names.push_back("grad(" + op_name + ")");
      prisms::ConditionalOStreams::pout_base() << op_name << " ";
    }
  prisms::ConditionalOStreams::pout_base() << "\n";

  // Assign fields
  int var_index = 0;
  for (const auto &mu_name : mu_names)
    {
      set_variable_name(var_index, mu_name);
      set_variable_type(var_index, Scalar);
      set_variable_equation_type(var_index, ExplicitTimeDependent);
      insert_dependencies_value_term_rhs(var_index, mu_names);
      insert_dependencies_value_term_rhs(var_index, grad_mu_names);
      insert_dependencies_gradient_term_rhs(var_index, mu_names);
      insert_dependencies_gradient_term_rhs(var_index, grad_mu_names);
      var_index++;
    }
  for (const auto &op_name : op_names)
    {
      set_variable_name(var_index, op_name);
      set_variable_type(var_index, Scalar);
      set_variable_equation_type(var_index, ExplicitTimeDependent);
      insert_dependencies_value_term_rhs(var_index, mu_names);
      insert_dependencies_value_term_rhs(var_index, grad_mu_names);
      insert_dependencies_value_term_rhs(var_index, op_names);
      insert_dependencies_value_term_rhs(var_index, grad_op_names);
      insert_dependencies_gradient_term_rhs(var_index, mu_names);
      insert_dependencies_gradient_term_rhs(var_index, grad_mu_names);
      insert_dependencies_gradient_term_rhs(var_index, op_names);
      insert_dependencies_gradient_term_rhs(var_index, grad_op_names);
      var_index++;
    }

  // Get names for comp fields
  std::vector<std::string> c_names;
  for (const auto &comp_name : system.comp_names)
    {
      c_names.push_back("c_" + comp_name);
    }
  // Assign fields
  for (const auto &c_name : c_names)
    {
      set_variable_name(var_index, c_name);
      set_variable_type(var_index, Scalar);
      set_variable_equation_type(var_index, ExplicitTimeDependent);
      set_is_postprocessed_field(var_index, true);
      insert_dependencies_value_term_rhs(var_index, mu_names);
      insert_dependencies_value_term_rhs(var_index, op_names);
      var_index++;
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block) const
{
  SystemContainer<dim, degree, number> sys(system_nd, get_user_inputs());
  uint                                 var_index = 0;
  sys.initialize_fields_explicit(variable_list, var_index);
  sys.calculate_locals();
  sys.calculate_detadt();
  sys.calculate_dmudt();
  var_index = 0;
  sys.submit_fields(variable_list, var_index);
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block,
  [[maybe_unused]] Types::Index                           index) const
{}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_lhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block,
  [[maybe_unused]] Types::Index                           index) const
{}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_postprocess_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block) const
{
  SystemContainer<dim, degree, number> ppsys(system, get_user_inputs());
  uint                                 var_index = 0;
  ppsys.initialize_fields_postprocess(variable_list, var_index);
  ppsys.calculate_sum_sq_eta();
  ppsys.calculate_h();
  uint pp_index = 0;
  ppsys.submit_pp_fields(variable_list, pp_index);
}

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE
