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


#ifndef RESIDUALVAR_H
#define RESIDUALVAR_H

#include "NodalKernel.h"
#include <petscsnes.h>

//Forward Declarations
class ResidualVar;

template<>
InputParameters validParams<ResidualVar>();

/**
 * Represents the rate in a simple ODE of du/dt = f
 */
class ResidualVar : public NodalKernel
{
public:
  /**
   * Constructor grabs the Function
   */
  ResidualVar(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
   const VariableValue & _res;
};

#endif
