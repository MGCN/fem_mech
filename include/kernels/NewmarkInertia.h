/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
/****************************************************************/	
/*               DO NOT MODIFY THIS HEADER                      */	
/*       MECH - ICS Mechanical simulation framework             */	
/*                Prepared by Maria Nestola,                    */	
/*                  ICS, USI, 6900 Lugano                       */	
/*                                                              */	
/*      Kernel to consider the intertial term.                  */
/*      It implements the term (_rho_s - _rho_f) * u_tt         */
/*      by using the Newmark method                             */
/****************************************************************/

#ifndef NewmarkInertia_H
#define NewmarkInertia_H

#include "Kernel.h"

//Forward Declarations
class NewmarkInertia;

template<>
InputParameters validParams<NewmarkInertia>();

class NewmarkInertia: public Kernel
{
public:

  NewmarkInertia(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();
    Real _rho_f;
    Real _rho_s;
    const VariableValue & _u_old;
    const VariableValue & _u_older;


};

#endif //NewmarkInertialH
