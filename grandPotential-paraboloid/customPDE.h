#include "ParaboloidSystem.h"
#include "SystemContainer.h"

#include <core/matrixFreePDE.h>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <string>

using namespace dealii;

template <int dim, int degree>
class customPDE : public MatrixFreePDE<dim, degree>
{
public:
  // Constructor
  customPDE(userInputParameters<dim> _userInputs)
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

  // Function to set the initial conditions (in ICs_and_BCs.h)
  void
  setInitialCondition([[maybe_unused]] const Point<dim>  &p,
                      [[maybe_unused]] const unsigned int index,
                      [[maybe_unused]] double            &scalar_IC,
                      [[maybe_unused]] Vector<double>    &vector_IC) override;

  // Function to set the non-uniform Dirichlet boundary conditions (in
  // ICs_and_BCs.h)
  void
  setNonUniformDirichletBCs([[maybe_unused]] const Point<dim>  &p,
                            [[maybe_unused]] const unsigned int index,
                            [[maybe_unused]] const unsigned int direction,
                            [[maybe_unused]] const double       time,
                            [[maybe_unused]] double            &scalar_BC,
                            [[maybe_unused]] Vector<double>    &vector_BC) override;

private:
#include <core/typeDefs.h>

  const userInputParameters<dim> userInputs;

  // Function to set the RHS of the governing equations for explicit time
  // dependent equations (in equations.h)
  void
  explicitEquationRHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;

  // Function to set the RHS of the governing equations for all other equations
  // (in equations.h)
  void
  nonExplicitEquationRHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;

  // Function to set the LHS of the governing equations (in equations.h)
  void
  equationLHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;

  // Function to set postprocessing expressions
  void
  postProcessedFields(
    [[maybe_unused]] const variableContainer<dim, degree, VectorizedArray<double>>
      &variable_list,
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                              &pp_variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;

// Function to set the nucleation probability (in nucleation.h)
#ifdef NUCLEATION_FILE_EXISTS
  double
  getNucleationProbability([[maybe_unused]] variableValueContainer variable_value,
                           [[maybe_unused]] double                 dV) const override;
#endif

  // ================================================================
  // Methods specific to this subclass
  // ================================================================
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
   * @brief Fraction of the theoretical maximum time step to use
   */
  double timestep_alpha = userInputs.get_model_constant_double("timestep_alpha");
  /**
   * @brief Map of the initial grand potential densities for each phase (used for
   * printing)
   */
  std::map<std::string, double> initial_omega;

  double r0 = userInputs.get_model_constant_double("r0");
  // ================================================================
};
