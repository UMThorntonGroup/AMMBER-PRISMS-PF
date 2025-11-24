// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <cmath>
#include <prismspf/config.h>
#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/nonuniform_dirichlet.h>
#include <prismspf/user_inputs/user_input_parameters.h>
#include <prismspf/utilities/utilities.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::set_initial_condition(
  [[maybe_unused]] const unsigned int       &index,
  [[maybe_unused]] const unsigned int       &component,
  [[maybe_unused]] const dealii::Point<dim> &point,
  [[maybe_unused]] number                   &scalar_value,
  [[maybe_unused]] number                   &vector_component_value) const
{
  std::map<uint, uint>        op_index;
  std::map<std::string, uint> comp_index;
  uint                        var_index = 0;
  for (const auto &comp_name : system_nd.comp_names)
    {
      comp_index[comp_name] = var_index++;
    }
  uint eta_index = 0;
  for ([[maybe_unused]] const auto &phase_name : system_nd.order_params)
    {
      op_index[eta_index++] = var_index++;
    }

  // ---------------------------------------------------------------------
  // TODO: ENTER THE INITIAL CONDITIONS HERE >
  // ---------------------------------------------------------------------
  // Custom coordinate system
  auto   domain_size = get_user_inputs().get_spatial_discretization().get_size();
  double center[3]   = {0.5 * domain_size[0],
                      (dim > 1) ? 0.5 * domain_size[1] : 0.0,
                      (dim > 2) ? 0.5 * domain_size[2] : 0.0};
  double x           = point[0] - center[0];
  double y           = (dim < 2) ? 0.0 : point[1] - center[1];
  double z           = (dim < 3) ? 0.0 : point[2] - center[2];
  double r2          = x * x + y * y + z * z;
  (void) r2; // To avoid unused variable warning

  // TODO: Make relevant geometries
  [[maybe_unused]] double circular = system_nd.interface(0.5 * (r0 * r0 - r2) / r0);
  [[maybe_unused]] double flat     = system_nd.interface(r0 - y);

  // TODO: Populate eta0 with the initial condition for the order parameters
  std::vector<double> eta0(system_nd.order_params.size(), 0.0);
  eta0[0] = 1.0 - circular;
  eta0[1] = circular;
  // ---------------------------------------------------------------------
  //  < ENTER THE INITIAL CONDITIONS HERE
  // ---------------------------------------------------------------------

  // Submit the fields
  var_index = 0;
  for (uint comp_index = 0; comp_index < system_nd.comp_names.size(); comp_index++)
    {
      if (index == var_index)
        {
          double mu0 = 0.;
          eta_index  = 0;
          for (const auto &phase_index : system_nd.order_params)
            {
              auto &phase_comp_info =
                system_nd.phases.at(phase_index).comps.at(comp_index);
              mu0 += eta0[eta_index] * phase_comp_info.k_well *
                     (phase_comp_info.x0 - phase_comp_info.c_min);
              eta_index++;
            }
          scalar_value = mu0;
          return;
        }
      var_index++;
    }
  eta_index = 0;
  for ([[maybe_unused]] const auto &phase_name : system_nd.order_params)
    {
      if (index == op_index[eta_index])
        {
          scalar_value = eta0[eta_index];
          return;
        }
      eta_index++;
      var_index++;
    }

  // ---------------------------------------------------------------------
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::set_nonuniform_dirichlet(
  [[maybe_unused]] const unsigned int       &index,
  [[maybe_unused]] const unsigned int       &boundary_id,
  [[maybe_unused]] const unsigned int       &component,
  [[maybe_unused]] const dealii::Point<dim> &point,
  [[maybe_unused]] number                   &scalar_value,
  [[maybe_unused]] number                   &vector_component_value) const
{
  this->set_initial_condition(index,
                              component,
                              point,
                              scalar_value,
                              vector_component_value);
}

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE