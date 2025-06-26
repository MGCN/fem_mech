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

#ifndef VALUEAUX_H
#define VALUEAUX_H

#include "AuxKernel.h"

//Forward Declarations
class ValueAux;


template<>
InputParameters validParams<ValueAux>();

/**
 * Function auxiliary value
 */
class ValueAux : public AuxKernel
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  ValueAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
  const VariableValue & _v;

  /// Function being used to compute the value of this kernel

};

#endif // FUNCTIONAUX_H
