// ===========================================================================
// FUNCTION FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::setInitialCondition([[maybe_unused]] const Point<dim>  &p,
                                            [[maybe_unused]] const unsigned int index,
                                            [[maybe_unused]] double            &scalar_IC,
                                            [[maybe_unused]] Vector<double>    &vector_IC)
{
  std::map<uint, uint>        op_var_index;
  std::map<std::string, uint> comp_index;
  uint                        var_index = 0;
  for (const auto &comp_name : isoSys.comp_names)
    {
      comp_index[comp_name] = var_index++;
    }
  uint eta_index = 0;
  for ([[maybe_unused]] const auto &phase_name : isoSys.order_params)
    {
      op_var_index[eta_index++] = var_index++;
    }

  // ---------------------------------------------------------------------
  // TODO: ENTER THE INITIAL CONDITIONS HERE >
  // ---------------------------------------------------------------------
  // Custom coordinate system
  double center[3] = {0.5 * userInputs.domain_size[0],
                      0.5 * userInputs.domain_size[1],
                      (dim > 2) * userInputs.domain_size[2]};
  double x         = p[0] - center[0];
  double y         = p[1] - center[1];
  double z         = (dim < 3) ? 0.0 : p[2] - center[2];
  double r2        = x * x + y * y + z * z;
  (void) r2;

  // TODO: Make relevant geometries
  [[maybe_unused]] double circular = interface(0.5 * (r0 * r0 - r2) / r0);
  [[maybe_unused]] double flat     = interface(0.5 * (r0 * r0 - y * y) / r0);

  // TODO: Populate eta0 with the initial condition for the order parameters
  std::vector<double> eta0(isoSys.order_params.size(), 0.0);
  eta0[0] = 1.0 - circular;
  eta0[1] = circular;
  // ---------------------------------------------------------------------
  // < ENTER THE INITIAL CONDITIONS HERE
  // ---------------------------------------------------------------------

  // Compute the initial condition for the chemical potential
  boost_vector<double> mu0(isoSys.num_comps);
  for (const auto &phase_index : isoSys.order_params)
    {
      auto &phase_info = isoSys.phases.at(phase_index);
      mu0 += phase_info.mu0 * eta0[phase_index];
    }
  // Submit the initial condition for the chemical potential
  var_index = 0;
  for (uint comp_index = 0; comp_index < isoSys.num_comps; comp_index++)
    {
      if (index == var_index)
        {
          scalar_IC = mu0[comp_index];
        }
      var_index++;
    }
  eta_index = 0;
  for ([[maybe_unused]] const auto &phase_name : isoSys.order_params)
    {
      if (index == op_var_index[eta_index])
        {
          scalar_IC = eta0[eta_index];
        }
      eta_index++;
      var_index++;
    }

  // ---------------------------------------------------------------------
}

// ===========================================================================
// FUNCTION FOR NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS
// ===========================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::setNonUniformDirichletBCs(
  [[maybe_unused]] const Point<dim>  &p,
  [[maybe_unused]] const unsigned int index,
  [[maybe_unused]] const unsigned int direction,
  [[maybe_unused]] const double       time,
  [[maybe_unused]] double            &scalar_BC,
  [[maybe_unused]] Vector<double>    &vector_BC)
{
  // --------------------------------------------------------------------------
  // ENTER THE NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS HERE
  // --------------------------------------------------------------------------
  // Use the initial condition function to set the boundary conditions
  this->setInitialCondition(p, index, scalar_BC, vector_BC);
  // --------------------------------------------------------------------------
}
