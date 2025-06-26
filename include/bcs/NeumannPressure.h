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

#ifndef NEUMANNPRESSURE_H
#define NEUMANNPRESSURE_H

#include "IntegratedBC.h"

// Forward Declarations
class NeumannPressure;
class Function;

template <>
InputParameters validParams<NeumannPressure>();

/**
 * Implements a simple constant Neumann BC where grad(u)=alpha * v on the
 * boundary.
 * Uses the term produced from integrating the diffusion operator by parts.
 */
class NeumannPressure : public IntegratedBC
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be
   * built using the same
   * constructor.
   */
  NeumannPressure(const InputParameters &parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  /**
   * Multiplier on the boundary.
   */
  Real _pressure;
  unsigned _component;
//  const VariableValue & _pressureFluid;
  const Function &_function2;

  const VariableGradient &_grad_disp_x;
  const VariableGradient &_grad_disp_y;
  const VariableGradient &_grad_disp_z;
};
#endif
