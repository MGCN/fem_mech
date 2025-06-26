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

#ifndef INCOMPRESSIBLESTRESSDIVERGENCELIBMESHTENSOR_H
#define INCOMPRESSIBLESTRESSDIVERGENCELIBMESHTENSOR_H

#include "Kernel.h"
#include "ElasticMaterial.h"

class IncompressibleStressDivergenceLibmeshTensor;

template <>
InputParameters validParams<IncompressibleStressDivergenceLibmeshTensor>();

class IncompressibleStressDivergenceLibmeshTensor : public Kernel
{
public:
  IncompressibleStressDivergenceLibmeshTensor(
      InputParameters const &parameters);

protected:
  virtual void computeJacobian();

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  unsigned int _component;
  unsigned const _numComp;

  unsigned const _p_var;

  RealVectorValue zeros;

  RealTensorValue H;

  const VariableValue &_p;

  MaterialProperty<Real> const &_J;
  MaterialProperty<RealTensorValue> const &_invFtr;

  const MaterialProperty<RealTensorValue> &_P;

  RealTensorValue **_stress_lin;

  const MaterialProperty<ElasticMaterial *> &_materialPointer;
};

#endif // INCOMPRESSIBLESTRESSDIVERGENCELIBMESHTENSOR_H
