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
/* Kernel to ensure conservation of momentum for compressible   */
/* hyperelastic materials.                                      */
/****************************************************************/

#ifndef STRESSDIVERGENCELIBMESHTENSOR_H
#define STRESSDIVERGENCELIBMESHTENSOR_H

#include "Kernel.h"
#include "ElasticMaterial.h"

class StressDivergenceLibmeshTensor;

template <>
InputParameters validParams<StressDivergenceLibmeshTensor>();

class StressDivergenceLibmeshTensor : public Kernel
{
public:
  StressDivergenceLibmeshTensor(InputParameters const &parameters);

  unsigned int _component;

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const MaterialProperty<RealTensorValue> &_P;

  RealTensorValue **_stress_lin;

  const MaterialProperty<ElasticMaterial *> &_materialPointer;

  unsigned const _numComp;

  RealVectorValue zeros;

  RealTensorValue H;
};

#endif // STRESSDIVERGENCELIBMESHTENSOR_H
