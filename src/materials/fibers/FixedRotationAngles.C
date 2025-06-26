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
/*                Prepared by Maria Nestola                     */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/* Material to compute fibers direction.                        */
/****************************************************************/


#include "FixedRotationAngles.h"

registerMooseObject("MECHApp", FixedRotationAngles);

template <>
InputParameters
validParams<FixedRotationAngles>()
{

   // inherit the parameters of the FibersMaterial:
  InputParameters params = validParams<FibersMaterial>();

  // specify sheet directions:
  params.addRequiredParam<Real>("beta0", "beta0");
  params.addRequiredParam<Real>("beta1", "beta1");

  return params;
}

FixedRotationAngles::FixedRotationAngles(const InputParameters &parameters):

   // inherit the parameters of the FibersMaterial:  
    FibersMaterial(parameters),
   
   // inherit the parameters of the FibersMaterial:
     beta0(getParam<Real>("beta0")),
     beta1(getParam<Real>("beta1"))
{




 // E_fiber0 = RealVectorValue(1, 0, 1);
  E_fiber = RealVectorValue(1.0 * std::cos(beta0 * M_PI / 180.0),   1.0 * std::sin(beta0 * M_PI / 180.0),0);
  E_sheet = RealVectorValue(1.0 * std::cos(beta1 * M_PI / 180.0), - 1.0 * std::sin(beta1 * M_PI / 180.0),0);
 // prima il secondo era -1
  
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
FixedRotationAngles::computeQpProperties()
{
  for (unsigned i = 0; i < 3; ++i)
  {
    _fiberDirection[_qp](i)  = E_fiber(i);
    _sheetDirection[_qp](i)  = E_sheet(i);
    _normalDirection[_qp](i) = E_normal(i);
  }

  computeQpTensorProperties();
}
