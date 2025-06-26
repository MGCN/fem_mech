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
/* Auxiliary Kernel to visualize the fibers.                    */
/****************************************************************/

#include "FibersAux.h"
registerMooseObject("MECHApp", FibersAux);


template <>
InputParameters
validParams<FibersAux>()
{
  // inherit the parameters of AuxKernel:
  InputParameters params = validParams<AuxKernel>();

  //specify for which component we consider:
  params.addRequiredParam<unsigned>("component", "component");
  
  return params;
}

FibersAux::FibersAux(const InputParameters &parameters): 

    // inherit the parameters of AuxKernel:
   AuxKernel(parameters),

  //specify for which component we consider:
   _component(getParam<unsigned>("component")),

   // inherit some material properties:
   _fiberDirection(getMaterialProperty<RealVectorValue>("fiberDirection"))
{
}

Real
FibersAux::computeValue()
{
//  std::cout<<"_fiberDirection[_qp](_component)==>"<<_fiberDirection[_qp](_component)<<std::endl;
  return _fiberDirection[_qp](_component);
}
