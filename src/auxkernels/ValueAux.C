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

#include "ValueAux.h"
registerMooseObject("MECHApp", ValueAux);

template<>
InputParameters validParams<ValueAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("v", "The variable whose value we are to match.");
  return params;
}

ValueAux::ValueAux(const InputParameters & parameters) :
    AuxKernel(parameters),
    _v(coupledValue("v"))
{
}

Real
ValueAux::computeValue()
{
    //std::cout<<"resn"<<_v[_qp]<<std::endl;
    return 1.0 * _v[_qp];
}

