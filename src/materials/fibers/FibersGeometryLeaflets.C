/****************************************************************/

/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*        MECH - ICS Imm Bound  simulation framework            */
/*                Prepared by Maria Nestola                      */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/*         Kernel for the Nearly Incompressibility              */
/****************************************************************/

#include "FibersGeometryLeaflets.h"

registerMooseObject("MECHApp", FibersGeometryLeaflets);

template <>
InputParameters
validParams<FibersGeometryLeaflets>()
{
  InputParameters params = validParams<FibersMaterial>();
  params.addRequiredCoupledVar(
      "thickness_parameter",
      "Aux variable with the thickness parameter (normalized "
      "mural position), will most often be produced with a "
      "CardiacThicknessParameterAux kernel. Make sure that "
      "its distinguishLVRV==false.");
  params.addRequiredParam<Real>("angle","angle between fiber families");
  return params;
}

FibersGeometryLeaflets::FibersGeometryLeaflets(const InputParameters &parameters)
    : FibersMaterial(parameters), _e(coupledValue("thickness_parameter")),
      _grad_e(coupledGradient("thickness_parameter")),
      _angle(getParam<Real>("angle"))
    
{
}

void
FibersGeometryLeaflets::computeQpProperties()
{

  
  const Real pi(3.141592653589);
  
  const Real R = _angle * pi / 180.0;

  const Real bracket(2. * _e[_qp] - 1.);
 
  Real alpha(R * bracket * bracket * bracket);


  // we already know the normal vector's direction (negative because of our e
  // being (1-e) of \ref Potse2006)
  RealVectorValue E_normal;
  
  if (_grad_e[_qp].norm() > 0)
  {
    E_normal = VectorNormalize(-_grad_e[_qp]);
  }
  else
  {
    // The gradient of the thickness parameter vanishes here.
    /// @todo TODO: can we find a better en in these cases than the spherical
    /// normal?
    E_normal = VectorNormalize(RealVectorValue(_q_point[_qp]));
  }

  const Real ez_en(E_normal(2)); ///< \f$\hat{e}_z\cdot\hat{e}_n\f$

  RealVectorValue ew;
  if (std::abs(ez_en) == 1.0)
  {
    // ez and en are (anti)parallel
    // for simplicity, we make ew the cylindrical normal vector here
    /// @todo TODO: can we find a better ew in these cases or average over
    /// neighbor cells, etc?
    ew = VectorNormalize(RealVectorValue(_q_point[_qp](0), _q_point[_qp](1), 0.));
  }
  else
  {
    const Real wz(1. / std::sqrt(1 - ez_en * ez_en));
    const Real wn(-wz * ez_en);
    ew = VectorNormalize(RealVectorValue(wn * E_normal(0), wn * E_normal(1),
                                         wn * E_normal(2) + wz));
  }
  // none of the VectorNormalize calls in the next lines should be necessary.
  // However, we want to avoid errors due to limited numeric precision.
  const RealVectorValue ev(VectorNormalize(ew.cross(E_normal)));
  
  E_fiber  = VectorNormalize(std::cos(R) * ev + std::sin(R) * ew);

  E_gfiber = VectorNormalize(std::cos(R) * ev - std::sin(R) * ew);
  
  E_sheet  = VectorNormalize(E_fiber.cross(E_normal));

  _normalDirection[_qp] = E_normal;

  _fiberDirection[_qp] = E_fiber;

  _gfiberDirection[_qp]  = E_gfiber;

  _sheetDirection[_qp] = E_sheet;

  computeQpTensorProperties();
}
