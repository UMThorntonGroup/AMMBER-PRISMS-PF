// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/parse_cmd_options.h>
#include <prismspf/core/problem.h>
#include <prismspf/utilities/logger.h>
#include <prismspf/utilities/utilities.h>

using namespace prismspf;

int
main(int argc, char *argv[])
{
  // Initialize MPI
  prismspf::MPIInitFinalize mpi_init(argc, argv);

  // Parse the command line options (if there are any) to get the name of the input
  // file
  ParseCMDOptions cli_options(argc, argv);

  constexpr unsigned int dim    = 2;
  constexpr unsigned int degree = 1;
  using number                  = float;

  ParaboloidSystem sys;
  std::ifstream    file("system.json");
  if (!file.is_open())
    throw std::runtime_error("Could not open system.json");
  nlohmann::json model_parameters;
  file >> model_parameters;
  sys.from_json(model_parameters);

  std::vector<FieldAttributes> field_attributes = sys.load_fields();
  std::vector<SolveBlock>      solve_blocks     = sys.load_blocks();

  // Set up user inputs
  UserInputParameters<dim>    user_inputs(cli_options.get_parameters_filename());
  SpatialDiscretization<dim> &space = user_inputs.spatial_discretization;
  TemporalDiscretization     &time  = user_inputs.temporal_discretization;

  // Choose the timestep automatically based on the CFL condition
  const double dx               = space.rectangular_mesh.size[0] / double(1 << space.global_refinement);
  const double stability_factor = user_inputs.user_constants.get_double("stability_factor");
  time.dt = stability_factor * prismspf::cfl_timestep<dim, degree>(sys.max_gradient_coefficient(), dx);
  Logger::instance() << "Set dt = " << time.dt << "\n";

  // Set up refinement criteria for the order parameters
  const RefinementCriterion op_criterion(RefinementFlags::Value, 0.01, 0.99);
  for (const auto &name : sys.get_order_parameter_names())
    {
      space.refinement_criteria["eta_" + name] = op_criterion;
    }
  Logger::instance() << "Set order parameter refinement criteria\n";

  // Set up the PDE operator and run the problem
  PhaseFieldTools<dim>           pf_tools;
  CustomPDE<dim, degree, number> pde_operator(user_inputs, pf_tools, sys);
  Problem<dim, degree, number>   problem(field_attributes, solve_blocks, user_inputs, pf_tools, pde_operator);
  problem.solve();

  return 0;
}
