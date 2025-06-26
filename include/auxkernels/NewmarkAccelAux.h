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
/*                                                              */		
/*          Auxkernel to compute acceleration field             */
/****************************************************************/

#ifndef NEWMARKACCELAUX_H
#define NEWMARKACCELAUX_H

#include "AuxKernel.h"

class NewmarkAccelAux;

template<>
InputParameters validParams<NewmarkAccelAux>();

class NewmarkAccelAux : public AuxKernel
{
public:

  /**
  *Computes Acceleration using Newmark Time integration scheme
  */
  NewmarkAccelAux(const InputParameters & parameters);

  virtual ~NewmarkAccelAux() {}

protected:
  virtual Real computeValue();

  const VariableValue & _disp_old;
  const VariableValue & _disp;
  const VariableValue & _vel_old;

  Real _beta;

  const VariableValue & u_old;

};

#endif //NEWMARKACCELAUX_H
