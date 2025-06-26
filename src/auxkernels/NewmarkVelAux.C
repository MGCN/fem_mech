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

#include "NewmarkVelAux.h"
registerMooseObject("MECHApp", NewmarkVelAux);

template<>
InputParameters validParams<NewmarkVelAux>()
{
    InputParameters params = validParams<AuxKernel>();
    params.addRequiredCoupledVar("acceleration","acceleration variable");
    params.addRequiredParam<Real>("gamma","gamma parameter");
    return params;
}

NewmarkVelAux::NewmarkVelAux(const InputParameters & parameters) :
AuxKernel(parameters),
_accel_old(coupledValueOld("acceleration")),
_accel(coupledValue("acceleration")),
_gamma(getParam<Real>("gamma")),
u_old(uOld())
{
}

Real
NewmarkVelAux::computeValue()
{
    Real vel_old = u_old[_qp];
    if (!isNodal())
        mooseError("must run on a nodal variable");
    return vel_old + ( _dt * ( 1 -_gamma)) * _accel_old[_qp] + _gamma * _dt * _accel[_qp];
}
