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
/*                Prepared by Marco Favino,                     */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/* Kernel to ensure the mass balance. To be used only coupled   */
/* with IncompressibleStressDivergenceLibmeshTensor             */
/****************************************************************/

#include "Mass.h"
#include "MooseMesh.h"

#include "libmesh/quadrature.h"
registerMooseObject("MECHApp", Mass);
template <>
InputParameters
validParams<Mass>()
{

  // inherit the parameters of the Kernels:
  InputParameters params = validParams<Kernel>();

  params.addParam<Real>("diffusionCoefficient", 1,"Reaction Coefficient");

  return params;
}

Mass::Mass(InputParameters const &parameters)
    :

      // inherit the parameters of the Kernels:
      Kernel(parameters),     
      _diffusionCoefficient(getParam<Real>("diffusionCoefficient"))

{

}

Real
Mass::computeQpResidual()
{
       return _diffusionCoefficient *_u[_qp] * _test[_i][_qp];

}

Real
Mass::computeQpJacobian()
{
  // the derivative of (J-1) w.r.t p
   return _diffusionCoefficient *_test[_i][_qp] * _phi[_j][_qp];
}

Real
Mass::computeQpOffDiagJacobian(unsigned int jvar)
{
      
    return 0.0;
}
