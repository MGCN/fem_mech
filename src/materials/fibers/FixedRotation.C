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

/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*       MECH - ICS Mechanical simulation framework             */
/*                Prepared by Marco Favino,                     */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/* Material to compute fibers direction.                        */
/****************************************************************/

#include "FixedRotation.h"

registerMooseObject("MECHApp", FixedRotation);

template <>
InputParameters
validParams<FixedRotation>()
{

   // inherit the parameters of the FibersMaterial:
  InputParameters params = validParams<FibersMaterial>();

  // specify sheet directions:
  params.addRequiredParam<RealVectorValue>("E_fiber", "E_fiber");
  params.addRequiredParam<RealVectorValue>("E_sheet", "E_sheet");

  return params;
}

FixedRotation::FixedRotation(const InputParameters &parameters):

   // inherit the parameters of the FibersMaterial:  
    FibersMaterial(parameters),
   
   // inherit the parameters of the FibersMaterial:
     E_fiber(getParam<RealVectorValue>("E_fiber")), E_sheet(getParam<RealVectorValue>("E_sheet"))
{

  //compute the normalized vectors
  Real norm_E_fiber =
      std::sqrt(E_fiber(0) * E_fiber(0) + E_fiber(1) * E_fiber(1) +
                E_fiber(2) * E_fiber(2));

  for (int i = 0; i < 3; ++i)
    E_fiber(i) = E_fiber(i) / norm_E_fiber;

  E_sheet = E_sheet - E_fiber.contract(E_sheet) * E_fiber;

  Real norm_E_sheet =
      std::sqrt(E_sheet(0) * E_sheet(0) + E_sheet(1) * E_sheet(1) +
                E_sheet(2) * E_sheet(2));

  for (int i = 0; i < 3; ++i)
    E_sheet(i) = E_sheet(i) / norm_E_sheet;

//compute the third remaining ortonormal vector 
  E_normal = E_fiber.cross(E_sheet);
}

void
FixedRotation::initQpStatefulProperties()
{

  //compute the normalized vectors
  Real norm_E_fiber =
      std::sqrt(E_fiber(0) * E_fiber(0) + E_fiber(1) * E_fiber(1) +
                E_fiber(2) * E_fiber(2));

  for (int i = 0; i < 3; ++i)
    E_fiber(i) = E_fiber(i) / norm_E_fiber;

  E_sheet = E_sheet - E_fiber.contract(E_sheet) * E_fiber;

  Real norm_E_sheet =
      std::sqrt(E_sheet(0) * E_sheet(0) + E_sheet(1) * E_sheet(1) +
                E_sheet(2) * E_sheet(2));

  for (int i = 0; i < 3; ++i)
    E_sheet(i) = E_sheet(i) / norm_E_sheet;

//compute the third remaining ortonormal vector 
  E_normal = E_fiber.cross(E_sheet);
}


void
FixedRotation::computeQpProperties()
{


  for (unsigned i = 0; i < 3; ++i)
  {
    _fiberDirection[_qp](i) = E_fiber(i);
    _sheetDirection[_qp](i) = E_sheet(i);
    _normalDirection[_qp](i) = E_normal(i);
  }





  computeQpTensorProperties();
}
