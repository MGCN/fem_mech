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

#ifndef FIXEDROTATION_H
#define FIXEDROTATION_H

#include "FibersMaterial.h"

class FixedRotation;

template <>
InputParameters validParams<FixedRotation>();

class FixedRotation : public FibersMaterial
{
public:
  FixedRotation(const InputParameters &parameters);

  void initQpStatefulProperties();
  void computeQpProperties();

  RealVectorValue E_fiber;
  RealVectorValue E_sheet;
  RealVectorValue E_normal;

};

#endif //FIXEDROTATION_H
