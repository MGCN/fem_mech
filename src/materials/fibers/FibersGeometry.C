#include "FibersGeometry.h"

registerMooseObject("MECHApp", FibersGeometry);

template <>
InputParameters
validParams<FibersGeometry>()
{
  InputParameters params = validParams<FibersMaterial>();
  params.addRequiredCoupledVar(
      "thickness_parameter",
      "Aux variable with the thickness parameter (normalized "
      "mural position), will most often be produced with a "
      "CardiacThicknessParameterAux kernel. Make sure that "
      "its distinguishLVRV==false.");
  return params;
}

FibersGeometry::FibersGeometry(const InputParameters &parameters)
    : FibersMaterial(parameters), _e(coupledValue("thickness_parameter")),
      _grad_e(coupledGradient("thickness_parameter"))
{
}

void
FibersGeometry::computeQpProperties()
{

  // initial angle for fibre angle as stated by \ref Potse2006
  // We follow the notation from this paper. Our e is defined slightly
  // different, though.
  const Real pi(3.141592653589);
  Real R =  1.0 * pi / 3.0;
  Real bracket = (2.0 * _e[_qp] - 1.0);

  if (_e[_qp] < 0.0)
  {
    R = pi / 4.0;
    bracket = -2.0 * _e[_qp] - 1.0;
  }


  Real alpha(R * bracket * bracket * bracket);

  // we already know the normal vector's direction (negative because of our e
  // being (1-e) of \ref Potse2006)
  RealVectorValue E_normal;
  // if (_grad_e[_qp].size() > 0)
  // {


     E_normal = VectorNormalize(_grad_e[_qp]);

    if (_e[_qp] > 0.0)
    {

      E_normal= - 1.0 * E_normal;

    }    


  // }
  // else
  // {
  //   // The gradient of the thickness parameter vanishes here.
  //   /// @todo TODO: can we find a better en in these cases than the spherical
  //   /// normal?
  //   E_normal = VectorNormalize(RealVectorValue(_q_point[_qp]));
  // }

  const Real ez_en(E_normal(2)); ///< \f$\hat{e}_z\cdot\hat{e}_n\f$

  RealVectorValue ew;
  if (std::abs(ez_en) == 1.0)
  {
    // ez and en are (anti)parallel
    // for simplicity, we make ew the cylindrical normal vector here
    /// @todo TODO: can we find a better ew in these cases or average over
    /// neighbor cells, etc?
    ew = VectorNormalize(
        RealVectorValue(_q_point[_qp](0), _q_point[_qp](1), 0.));
  }
  else
  {
    const Real wz(-1. / std::sqrt(1.0 - ez_en * ez_en)); // before it was +1
    const Real wn(-wz * ez_en);
    ew = VectorNormalize(RealVectorValue(wn * E_normal(0), wn * E_normal(1),
                                         wn * E_normal(2) + wz));
  }
  // none of the VectorNormalize calls in the next lines should be necessary.
  // However, we want to avoid errors due to limited numeric precision.
  const RealVectorValue ev(VectorNormalize(ew.cross(E_normal)));
  E_fiber = VectorNormalize(std::cos(alpha) * ev + std::sin(alpha) * ew);
  E_sheet = VectorNormalize(E_fiber.cross(E_normal));

  _fiberDirection[_qp]  = E_fiber;
  _normalDirection[_qp] = E_normal;
  _sheetDirection[_qp]  = E_sheet;

  if (isnan(_fiberDirection[_qp](0)) || isnan(_fiberDirection[_qp](1)) || isnan(_fiberDirection[_qp](2)) )
  {
    std::cout<<"nan _fiberDirection ";
    std::cout<<_fiberDirection[_qp]<<std::endl;


  }


  if (isnan(_normalDirection[_qp](0)) || isnan(_normalDirection[_qp](1)) || isnan(_normalDirection[_qp](2)) )
  {
    std::cout<<"nan _normalDirection ";
    std::cout<<_normalDirection[_qp]<<std::endl;
    std::cout<<"position "<<_q_point[_qp]<<std::endl;

  }


  if (isnan(_sheetDirection[_qp](0)) || isnan(_sheetDirection[_qp](1)) || isnan(_sheetDirection[_qp](2)) )
  {
    std::cout<<"nan _sheetDirection ";
    std::cout<<_sheetDirection[_qp]<<std::endl;
  }






  computeQpTensorProperties();
}

