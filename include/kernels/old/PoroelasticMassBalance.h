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
/* Kernel to ensure poroelastic mass balance.                   */
/****************************************************************/

#ifndef POROELASTICMASSBALANCE_H
#define POROELASTICMASSBALANCE_H

#include "TimeKernel.h"
#include "ElasticMaterial.h"

class PoroelasticMassBalance;

template <>
InputParameters validParams<PoroelasticMassBalance>();

class PoroelasticMassBalance : public TimeKernel
{
public:
  PoroelasticMassBalance(InputParameters const &parameters);

  unsigned int _component;

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  unsigned const _numComp;

  RealVectorValue zeros;

  RealTensorValue H;

  MaterialProperty<Real> const &_J;
  MaterialProperty<Real> const &_J_old;

  MaterialProperty<RealTensorValue> const &_invFtr;
};

#endif // POROELASTICMASSBALANCE_H
