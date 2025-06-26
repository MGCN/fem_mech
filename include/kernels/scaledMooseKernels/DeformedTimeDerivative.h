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

#ifndef DEFORMEDTIMEDERIVATIVE_H
#define DEFORMEDTIMEDERIVATIVE_H

#include "TimeDerivative.h"
#include "Kernel.h"
#include "ElasticMaterial.h"

class DeformedTimeDerivative;

template <>
InputParameters validParams<DeformedTimeDerivative>();

class DeformedTimeDerivative : public TimeDerivative
{
public:
  DeformedTimeDerivative(const InputParameters &parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  unsigned const _numComp;

  RealVectorValue zeros;

  RealTensorValue H;

  const MaterialProperty<Real> &_J;
  MaterialProperty<RealTensorValue> const &_invFtr;

 double _surface_to_volume;
  double _capacitance;

  Real **_J_lin;
};

#endif
