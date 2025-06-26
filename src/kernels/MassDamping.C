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

#include "MassDamping.h"
#include "MooseMesh.h"

#include "libmesh/quadrature.h"

registerMooseObject("MECHApp", MassDamping);

template <>
InputParameters
validParams<MassDamping>()
{

  // inherit the parameters of the Kernels:
  InputParameters params = validParams<Kernel>();


  params.addRequiredParam<Real>("Damping", "Damping");

  return params;
}

MassDamping::MassDamping(InputParameters const &parameters)
    :

      // inherit the parameters of the Kernels:
      Kernel(parameters),
      _u_old(valueOld()),
      _D(getParam<Real>("Damping"))

{

}

Real
MassDamping::computeQpResidual()
{
       return _D * ( _u[_qp]  - _u_old[_qp]) /_dt * _test[_i][_qp];

}

Real
MassDamping::computeQpJacobian()
{
    return _D * _test[_i][_qp] * _phi[_j][_qp] / _dt;
}

Real
MassDamping::computeQpOffDiagJacobian(unsigned int jvar)
{
      
    return 0.0;
}
