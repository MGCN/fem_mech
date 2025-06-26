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
/*                                                              */	
/*      Kernel to consider the intertial term.                  */
/*      It implements the term                                  */
/*     _int ((_rho_s - _rho_f) * u_tt * test * dV )             */
/*      by using the Newmark method                             */
/****************************************************************/

#include "SplitNewmarkInertia.h"
#include "SubProblem.h"

registerMooseObject("MECHApp", SplitNewmarkInertia);

template<>
InputParameters validParams<SplitNewmarkInertia>()
{
  InputParameters params = validParams<Kernel>();
  // params.addRequiredParam<Real>("rho_f", "density_fluid");
  // params.addRequiredParam<Real>("rho_s", "density_solid");
  params.addRequiredCoupledVar("velocity","velocity variable");
  params.addRequiredCoupledVar("rho","density");

  return params;
}

SplitNewmarkInertia::SplitNewmarkInertia(const InputParameters & parameters) :
    Kernel(parameters),
    // _rho_f(getParam<Real>("rho_f")),
    // _rho_s(getParam<Real>("rho_s")),
    _rho(coupledValue("rho")),
    _vel(coupledValueOld("velocity")),
    _u_old(valueOld()),
    _u_older(valueOlder())

{}

Real
SplitNewmarkInertia::computeQpResidual()
{
    // Real diff_rho = _rho_s - _rho_f;
    

    return _rho[_qp] * 4 * (( _u[_qp] - _u_old[_qp])/_dt/_dt - _vel[_qp]/_dt) * _test[_i][_qp];
    
}

Real
SplitNewmarkInertia::computeQpJacobian()
{
    // Real diff_rho =_rho_s - _rho_f;

    return 0.0;
    //diff_rho * 4 * _phi[_j][_qp]/_dt/_dt * _test[_i][_qp];

}
