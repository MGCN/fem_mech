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

#include "GFibersAux.h"
registerMooseObject("MECHApp", GFibersAux);


template <>
InputParameters
validParams<GFibersAux>()
{
  // inherit the parameters of AuxKernel:
  InputParameters params = validParams<AuxKernel>();

  //specify for which component we consider:
  params.addRequiredParam<unsigned>("component", "component");
  
  return params;
}

GFibersAux::GFibersAux(const InputParameters &parameters): 

    // inherit the parameters of AuxKernel:
   AuxKernel(parameters),

  //specify for which component we consider:
   _component(getParam<unsigned>("component")),

   // inherit some material properties:
   _gfiberDirection(getMaterialProperty<RealVectorValue>("gfiberDirection"))
{
}

Real
GFibersAux::computeValue()
{
//  std::cout<<"and _gfiberDirection[_qp](_component)==>"<<_gfiberDirection[_qp](_component)<<std::endl;
  return _gfiberDirection[_qp](_component);
}
