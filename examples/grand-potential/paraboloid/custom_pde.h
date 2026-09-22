#include <ammber/core/paraboloid_pde.h>
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
  using GrandPotentialPDE<dim, degree, number>::submit_ic_from_fields;
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
    // ---------------------------------------------------------------------
    // TODO: ENTER THE INITIAL CONDITIONS HERE >
    // ---------------------------------------------------------------------
    // Custom coordinate system
    const dealii::Tensor<1, dim> &mesh_size = get_user_inputs().spatial_discretization.rectangular_mesh.size;
    double                        x         = (dim < 1) ? 0.0 : point[0];
    double                        y         = (dim < 2) ? 0.0 : point[1];
    double                        z         = (dim < 3) ? 0.0 : point[2];
    double                        r2        = x * x + y * y + z * z;

    // TODO: Make relevant geometries
    [[maybe_unused]] double circular = sys.interface(0.5 * (r0 * r0 - r2) / r0);
    [[maybe_unused]] double flat     = sys.interface(0.5 * (r0 * r0 - y * y) / r0);

    // TODO: Populate eta0 with the initial condition for the order parameters
    std::vector<number> eta0(sys.num_ops(), 0.0);
    eta0[0] = 1.0 - circular;
    eta0[1] = circular;

    // Submit eta0 to utility that sets ICs
    submit_ic_from_fields(eta0, index, scalar_value);
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
};

PRISMS_PF_END_NAMESPACE
