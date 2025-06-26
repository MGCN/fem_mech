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
/*      Kernel to consider the intertial term.                  */
/*      It implements the term (_rho_s - _rho_f) * u_tt         */
/*      by using the Newmark method                             */
/****************************************************************/

#ifndef SplitNewmarkInertia_H
#define SplitNewmarkInertia_H

#include "Kernel.h"

//Forward Declarations
class SplitNewmarkInertia;

template<>
InputParameters validParams<SplitNewmarkInertia>();

class SplitNewmarkInertia: public Kernel
{
public:

  SplitNewmarkInertia(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();
    // Real _rho_f;
    // Real _rho_s;
    const VariableValue & _rho;
    const VariableValue & _vel;
    const VariableValue & _u_old;
    const VariableValue & _u_older;
    




};

#endif //TwoStepsNewmarkInertialH
