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

#ifndef GROWTHFACTOR_H
#define GROWTHFACTOR_H

#include "Kernel.h"


class GrowthFactor;

template <>
InputParameters validParams<GrowthFactor>();

class GrowthFactor: public Kernel
{
public:
  GrowthFactor(InputParameters const &parameters);

    
  //const VariableValue & _u_old;
  Real _k_plus;
  Real _k_minus;
  
  const VariableValue & _threshold;

  Real _theta_plus;
  Real _theta_minus;
  Real _m_plus;
  Real _m_minus;

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);


};

#endif // MASS_H
