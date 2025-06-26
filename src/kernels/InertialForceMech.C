//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html


#include "InertialForceMech.h"

#include "SubProblem.h"

registerMooseObject("MECHApp", InertialForceMech);

template <>
InputParameters validParams<InertialForceMech>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Calculates the residual for the interial force " "($M \\cdot acceleration$) and the contribution of mass"
                             " dependent Rayleigh damping" " integration scheme ($\\eta \\cdot M \\cdot velq2");
  params.addParam<Real>( "density", 0., "Name of the density");
  params.addParam<Real>("eta", 0., "Rayleigh damping parameter");
  return params;
}

InertialForceMech::InertialForceMech(const InputParameters & parameters)
: Kernel(parameters),
_density(getParam<Real>("density")),
_u_old(valueOld()),
_u_older(valueOlder()),
_eta(getParam<Real>("eta"))
{
}
	
Real
InertialForceMech::computeQpResidual()
{
  if (_dt == 0)
      return 0;
  else
  {
    Real accel = 2.0*_u[_qp]/(_dt*_dt + _dt*_dt_old) - 2.0*_u_old[_qp]/(_dt*_dt_old) + 2.0*_u_older[_qp]/(_dt*_dt_old + _dt_old*_dt_old);
    // Real vel = (_u[_qp] - _u_old[_qp])/_dt;
    
    Real vel = _u[_qp]*(2.*_dt + _dt_old)/(_dt*_dt + _dt*_dt_old) - _u_old[_qp]*(_dt + _dt_old)/(_dt*_dt_old) + _u_older[_qp]*_dt/(_dt*_dt_old + _dt_old*_dt_old);
   
    return _test[_i][_qp] * _density * (accel + _eta * vel);
  }

}
	
Real
InertialForceMech::computeQpJacobian()
{
  if (_dt == 0)
  return 0;
  else
  return _test[_i][_qp] * _density * 2.0 / (_dt * _dt + _dt * _dt_old) * _phi[_j][_qp] + _test[_i][_qp] * _eta * _density * (2.*_dt + _dt_old)/(_dt*_dt + _dt*_dt_old) * _phi[_j][_qp];
// + _test[_i][_qp] * _eta * _density / _dt * _phi[_j][_qp];
}


