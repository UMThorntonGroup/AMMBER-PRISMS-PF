#include "paraboloid_pde.h"

#include <prismspf/core/pde_operator_base.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class CustomPDE : public GrandPotentialPDE<dim, degree, number>
{
public:
  using ScalarValue = dealii::VectorizedArray<number>;
  using ScalarGrad  = dealii::Tensor<1, dim, ScalarValue>;
  using ScalarHess  = dealii::Tensor<2, dim, ScalarValue>;
  using VectorValue = dealii::Tensor<1, dim, ScalarValue>;
  using VectorGrad  = dealii::Tensor<2, dim, ScalarValue>;
  using VectorHess  = dealii::Tensor<3, dim, ScalarValue>;
  using GrandPotentialPDE<dim, degree, number>::sys;
  using PDEOperatorBase<dim, degree, number>::get_user_inputs;
  using PDEOperatorBase<dim, degree, number>::get_pf_tools;

  double r0;

  /**
   * @brief Constructor.
   */
  CustomPDE(const UserInputParameters<dim> &_user_inputs,
            PhaseFieldTools<dim>           &_pf_tools,
            const ParaboloidSystem         &_sys)
    : GrandPotentialPDE<dim, degree, number>(_user_inputs, _pf_tools, _sys)
    , r0(get_user_inputs().user_constants.get_double("r0"))
  {}

  /**
   * @brief User-implemented class for the setting initial conditions.
   */
  void
  set_initial_condition([[maybe_unused]] const unsigned int       &index,
                        [[maybe_unused]] const unsigned int       &component,
                        [[maybe_unused]] const dealii::Point<dim> &point,
                        [[maybe_unused]] number                   &scalar_value,
                        [[maybe_unused]] number                   &vector_component_value) const override
  {
    const dealii::Tensor<1, dim> &mesh_size = get_user_inputs().spatial_discretization.rectangular_mesh.size;

    std::map<uint, uint>        op_index;
    std::map<std::string, uint> comp_index;
    uint                        var_index = 0;
    for (const auto &comp_name : sys.comp_names)
      {
        comp_index[comp_name] = var_index++;
      }
    uint eta_index = 0;
    for ([[maybe_unused]] const auto &phase_name : sys.order_params)
      {
        op_index[eta_index++] = var_index++;
      }

    // ---------------------------------------------------------------------
    // TODO: ENTER THE INITIAL CONDITIONS HERE >
    // ---------------------------------------------------------------------
    // Custom coordinate system
    const dealii::Point<dim> center(mesh_size * 0.5);
    double                   x  = (dim < 1) ? 0.0 : point[0] - center[0];
    double                   y  = (dim < 2) ? 0.0 : point[1] - center[1];
    double                   z  = (dim < 3) ? 0.0 : point[2] - center[2];
    double                   r2 = x * x + y * y + z * z;
    (void) r2;

    // TODO: Make relevant geometries
    [[maybe_unused]] double circular = sys.interface(0.5 * (r0 * r0 - r2) / r0);
    [[maybe_unused]] double flat     = sys.interface(0.5 * (r0 * r0 - y * y) / r0);

    // TODO: Populate eta0 with the initial condition for the order parameters
    std::vector<double> eta0(sys.num_ops(), 0.0);
    eta0[0] = 1.0 - circular;
    eta0[1] = circular;
    // ---------------------------------------------------------------------
    //  < ENTER THE INITIAL CONDITIONS HERE
    // ---------------------------------------------------------------------

    // Submit the fields
    const double sum_sq_eta = sum_sq(eta0) + 1e-8;
    for (uint comp_index = 0; comp_index < sys.num_comps(); comp_index++)
      {
        var_index = sys.mu_base() + comp_index;
        if (index == var_index)
          {
            double mu0 = 0.;
            eta_index  = 0;
            for (const auto &phase_index : sys.order_params)
              {
                auto &phase_comp_info = sys.phases.at(phase_index).comps.at(comp_index);
                mu0 += eta0[eta_index] * eta0[eta_index] / sum_sq_eta * phase_comp_info.k_well *
                       (phase_comp_info.x0 - phase_comp_info.c_min);
                eta_index++;
              }
            scalar_value = mu0;
            return;
          }
      }
    var_index = sys.eta_base();
    eta_index = 0;
    for ([[maybe_unused]] const auto &phase_name : sys.order_params)
      {
        if (index == op_index[eta_index])
          {
            scalar_value = eta0[eta_index];
            return;
          }
        eta_index++;
        var_index++;
      }
  }

  /**
   * @brief User-implemented function for setting Dirichlet boundary conditions. Default
   * behavior is to call initial conditions.
   */
  /* virtual void
  set_dirichlet([[maybe_unused]] const unsigned int       &index,
                [[maybe_unused]] const unsigned int       &boundary_id,
                [[maybe_unused]] const unsigned int       &component,
                [[maybe_unused]] const dealii::Point<dim> &point,
                [[maybe_unused]] const SimulationTimer    &sim_timer,
                [[maybe_unused]] number                   &scalar_value,
                [[maybe_unused]] number                   &vector_component_value) const override; */

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
