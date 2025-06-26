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
/*       MECH - ICS Mechanical simulation framework             */
/*                Prepared by Maria Nestola,                    */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/****************************************************************/

#include "ElementIntegralEnergyNeoHookean.h"
#include "Material.h"
#include "MooseMesh.h"
#include <math.h>

registerMooseObject("MECHApp", ElementIntegralEnergyNeoHookean);

template <>
InputParameters
validParams<ElementIntegralEnergyNeoHookean>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredParam<Real>("mu", "mu");
  params.addRequiredParam<Real>("kappa", "kappa");
  params.addRequiredCoupledVar("pressure", "pressure");
  return params;
}

ElementIntegralEnergyNeoHookean::ElementIntegralEnergyNeoHookean(
    const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    _mu(getParam<Real>("mu")),
    _kappa(getParam<Real>("kappa")),
    _J(getMaterialProperty<Real>("deformationDeterminant")),
    _C(getMaterialProperty<RealTensorValue>("CauchyGreen")),
    _J23(getMaterialProperty<Real>("J23")),
    _F(getMaterialProperty<RealTensorValue>("deformationGradient")),
    _pressure(coupledValue("pressure"))

{

}

Real
ElementIntegralEnergyNeoHookean::computeQpIntegral()
{
    
    Real W1(0.0);
    W1 =  _mu * (_J23[_qp] * _C[_qp].tr() - _mesh.dimension()) + 0.5 * _kappa/4.0 * ( (_J[_qp] - 1 ) * (_J[_qp] - 1 ) +  std::log(_J[_qp]) * std::log(_J[_qp])) + _pressure[_qp] * (_J[_qp] - 1 ) ;
    
  //  Real W2(0.0);
  //  W2 =  _mu * (_C[_qp].tr() - _mesh.dimension() - std::log(_J[_qp])) + _kappa *  std::log(_J[_qp]) *  std::log(_J[_qp])  ;
    
  return W1;
 }
