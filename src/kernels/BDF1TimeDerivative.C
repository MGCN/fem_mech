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
/*                Prepared by Maria Nestola                     */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/*         Kernel for the Nearly Incompressibility              */
/****************************************************************/

#include "BDF1TimeDerivative.h"

registerMooseObject("MECHApp", BDF1TimeDerivative);

template<>
InputParameters validParams<BDF1TimeDerivative>()
{
    InputParameters params = validParams<TimeDerivative>();
    params.addRequiredParam<Real>("rho_f", "density_fluid");
    params.addRequiredParam<Real>("rho_s", "density_solid");
    
    return params;
}

BDF1TimeDerivative::BDF1TimeDerivative(const InputParameters & parameters) :
TimeDerivative(parameters),
_u_old(valueOld()),
_u_older(valueOlder()),
_dt(_fe_problem.dt()),
_t_step(_fe_problem.timeStep()),
_rho_s(getParam<Real>("rho_s")),
_rho_f(getParam<Real>("rho_f"))

{}

Real
BDF1TimeDerivative::computeQpResidual()
{

    if (_t_step == 1 )
       
       return  (_rho_s - _rho_f) * _test[_i][_qp] * ( 1. / _dt / _dt) * ( _u[_qp] - _u_old[_qp]);
    
    else
        
        return (_rho_s - _rho_f) * _test[_i][_qp] * ( 1. / _dt / _dt) * (1. * _u[_qp] - 2 * _u_old[_qp] + 1 * _u_older[_qp]);
   
}

Real
BDF1TimeDerivative::computeQpJacobian()
{
    if (_t_step == 1 )
         return (_rho_s - _rho_f) * _test[_i][_qp] * _phi[_j][_qp] * (1. /_dt / _dt);
    else
         return (_rho_s - _rho_f) * 1. * _test[_i][_qp] * _phi[_j][_qp] * (1. / _dt /_dt);
        
}
