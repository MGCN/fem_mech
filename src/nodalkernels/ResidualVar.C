/****************************************************************/
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

#include "ResidualVar.h"
#include "MooseMesh.h"

registerMooseObject("MECHApp", ResidualVar);

template<>
InputParameters validParams<ResidualVar>()
{
  InputParameters params = validParams<NodalKernel>();
  params.addRequiredCoupledVar("res","res");
  return params;
}

ResidualVar::ResidualVar(const InputParameters & parameters) :
    NodalKernel(parameters),
    _res(coupledValue("res"))

{
}

Real
ResidualVar::computeQpResidual()
{
   return - _res[0];
}
