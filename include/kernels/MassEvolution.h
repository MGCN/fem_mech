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
/****************************************************************/

#ifndef MASSEVOLUTION_H
#define MASSEVOLUTION_H

#include "Kernel.h"


class MassEvolution;

template <>
InputParameters validParams<MassEvolution>();

class MassEvolution: public Kernel
{
public:
  MassEvolution(InputParameters const &parameters);

  Real _rho_0;

  const VariableValue & _theta;

  const VariableValue & _theta_old;


protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  // virtual Real computeQpOffDiagJacobian(unsigned int jvar);


};

#endif // MASS_H
