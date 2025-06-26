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
/****************************************************************/

#ifndef DEFORMEDUSERFORCINGFUNCTION_H
#define DEFORMEDUSERFORCINGFUNCTION_H

#include "Kernel.h"

// Forward Declarations
class DeformedUserForcingFunction;

template <>
InputParameters validParams<DeformedUserForcingFunction>();

class DeformedUserForcingFunction : public Kernel
{
public:
  DeformedUserForcingFunction(const InputParameters &parameters);

protected:
  Real f();

  /**
   * Computes test function * forcing function.
   */
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  unsigned const _numComp;

  RealVectorValue zeros;

  RealTensorValue H;

  const Function &_func;

  const MaterialProperty<Real> &_J;
  MaterialProperty<RealTensorValue> const &_invFtr;

  Real **_J_lin;
};

#endif // DEFORMEDUSERFORCINGFUNCTION_H
