// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================

void
customAttributeLoader::loadVariableAttributes()
{
// Load the model parameters
#include "ParaboloidSystem.h"

#include <fstream>
#include <iostream>
  nlohmann::json   model_parameters;
  ParaboloidSystem isoSys;
  std::ifstream    system_file;
  system_file.open("system.json");
  system_file >> model_parameters;
  system_file.close();
  isoSys.from_json(model_parameters);
  uint var_index = 0;
  isoSys.print_parameters();
  isoSys.load_variables(this, var_index);
}

// =================================================================================
// Set the attributes of the postprocessed field variables
// =================================================================================

void
customAttributeLoader::loadPostProcessorVariableAttributes()
{
#include "ParaboloidSystem.h"

#include <fstream>
#include <iostream>
  nlohmann::json   model_parameters;
  ParaboloidSystem isoSys;
  std::ifstream    system_file;
  system_file.open("system.json");
  system_file >> model_parameters;
  system_file.close();
  isoSys.from_json(model_parameters);
  uint pp_index = 0;
  isoSys.load_pp_variables(this, pp_index);
}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time
// dependent)
// =============================================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::explicitEquationRHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
  SystemContainer<dim, degree> sys(isoSys, userInputs);
  uint                         var_index = 0;
  sys.initialize_fields_explicit(variable_list, var_index);

  sys.calculate_deltas();
  sys.calculate_sum_sq_eta();
  sys.calculate_h();
  sys.calculate_dhdeta();
  sys.calculate_local_mobility();

  sys.calculate_dmudt();
  var_index = 0;
  sys.submit_fields(variable_list, var_index);
}

// =============================================================================================
// postProcessedFields: Set the postprocessing expressions
// =============================================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::postProcessedFields(
  [[maybe_unused]] const variableContainer<dim, degree, VectorizedArray<double>>
    &variable_list,
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                            &pp_variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
  [[maybe_unused]] const VectorizedArray<double>             element_volume) const
{
  SystemContainer<dim, degree> ppsys(isoSys, userInputs);
  uint                         var_index = 0;
  ppsys.initialize_fields_postprocess(variable_list, var_index);
  ppsys.calculate_deltas();
  ppsys.calculate_sum_sq_eta();
  ppsys.calculate_h();
  uint pp_index = 0;
  ppsys.submit_pp_fields(pp_variable_list, pp_index);
}

// =============================================================================================
// nonExplicitEquationRHS (needed only if one or more equation is time
// independent or auxiliary)
// =============================================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::nonExplicitEquationRHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{
  SystemContainer<dim, degree> sys(isoSys, userInputs);
  uint                         var_index = 0;
  sys.initialize_fields_nonexplicit(variable_list, var_index);

  sys.calculate_deltas();
  sys.calculate_omega_phase();
  sys.calculate_sum_sq_eta();
  sys.calculate_h();
  sys.calculate_dhdeta();

  sys.calculate_detadt();
  var_index = 0;
  sys.submit_aux_fields(variable_list, var_index);
}

// =============================================================================================
// equationLHS (needed only if at least one equation is time independent)
// =============================================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::equationLHS(
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
  [[maybe_unused]] const VectorizedArray<double> element_volume) const
{}