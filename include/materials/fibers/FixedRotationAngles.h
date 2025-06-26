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


#ifndef FixedRotationAngles_H
#define FixedRotationAngles_H

#include "FibersMaterial.h"

class FixedRotationAngles;

template <>
InputParameters validParams<FixedRotationAngles>();

class FixedRotationAngles : public FibersMaterial
{
public:
  FixedRotationAngles(const InputParameters &parameters);

  void computeQpProperties();
  Real beta0;
  Real beta1;
  RealVectorValue E_fiber;
  RealVectorValue E_sheet;
  RealVectorValue E_normal;
};

#endif //FixedRotationHolzapfel2000_H
