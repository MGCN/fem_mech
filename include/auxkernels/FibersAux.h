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


#ifndef FIBERSAUX_H
#define FIBERSAUX_H

#include "AuxKernel.h"

class FibersAux;

template <>
InputParameters validParams<FibersAux>();

class FibersAux : public AuxKernel
{
public:
  FibersAux(const InputParameters &parameters);

protected:
  virtual Real computeValue();

  unsigned _component;

  MaterialProperty<RealVectorValue> const &_fiberDirection;
};
#endif
