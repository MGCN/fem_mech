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
/*                                                              */
/*            Auxkernel to compute velocity field               */
/****************************************************************/

#include "EulerVelAux.h"

registerMooseObject("MECHApp", EulerVelAux);

template<>
InputParameters validParams<EulerVelAux>()
{
    InputParameters params = validParams<AuxKernel>();
    params.addRequiredCoupledVar("disp","disp variable");
    params.addRequiredParam<Real>("gamma","gamma parameter");
    return params;
}
EulerVelAux::EulerVelAux(const InputParameters & parameters) :
AuxKernel(parameters),
_u_old(coupledValueOld("disp")),
_u(coupledValue("disp")),
_gamma(getParam<Real>("gamma"))
{
}

Real
EulerVelAux::computeValue()
{
    Real vel_old = _u_old[_qp];
    if (!isNodal())
        mooseError("must run on a nodal variable");
    return (_u[_qp] - _u_old[_qp])/_dt;
//vel_old + ( _dt * ( 1 -_gamma)) * _accel_old[_qp] + _gamma * _dt * _accel[_qp];
}
