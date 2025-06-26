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

#include "MassEvolution.h"
#include "SubProblem.h"
#include "MooseMesh.h"

registerMooseObject("MECHApp", MassEvolution);

template<>
InputParameters validParams<MassEvolution>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<Real>("rho_0", "density");
  params.addRequiredCoupledVar("theta","theta");


    

  return params;
}

MassEvolution::MassEvolution(const InputParameters & parameters) :
    Kernel(parameters),
    _rho_0(getParam<Real>("rho_0")),
    _theta(coupledValue("theta")),
    _theta_old(coupledValueOld("theta"))

{}

Real
MassEvolution::computeQpResidual()
{
   
   return 3.0 * _theta[_qp] * _theta[_qp] * _rho_0 * (_theta[_qp] - 1.0 * _theta_old[_qp])/_dt;
}

Real
MassEvolution::computeQpJacobian()
{
    
    return 0.0;

}
