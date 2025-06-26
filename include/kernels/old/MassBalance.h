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
/* Kernel to ensure the mass balance. To be used only coupled   */
/* with IncompressibleStressDivergenceLibmeshTensor             */
/****************************************************************/

#ifndef MASSBALANCE_H
#define MASSBALANCE_H

#include "Kernel.h"
#include "ElasticMaterial.h"

class MassBalance;

template <>
InputParameters validParams<MassBalance>();

class MassBalance : public Kernel
{
public:
  MassBalance(InputParameters const &parameters);

  unsigned int _component;

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  unsigned const _numComp;

  RealVectorValue zeros;

  RealTensorValue H;

  MaterialProperty<Real> const &_J;
  MaterialProperty<RealTensorValue> const &_invFtr;

  Real **_J_lin;
};

#endif // MASSBALANCE_H
