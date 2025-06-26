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

/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*       MECH - ICS Mechanical simulation framework             */
/*                Prepared by Marco Favino,                     */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/* Fibers Material to compute fibers direction.                 */
/****************************************************************/

#ifndef FIBERSMATERIAL_H
#define FIBERSMATERIAL_H

#include "Material.h"

class FibersMaterial;

template <>
InputParameters validParams<FibersMaterial>();

class FibersMaterial : public Material
{
public:
  FibersMaterial(const InputParameters &parameters);
  virtual void computeQpProperties() = 0;
  virtual void computeQpTensorProperties();

protected:

  MaterialProperty<RealVectorValue> &_fiberDirection;
  MaterialProperty<RealVectorValue> &_gfiberDirection;
  MaterialProperty<RealVectorValue> &_sheetDirection;
  MaterialProperty<RealVectorValue> &_normalDirection;
  MaterialProperty<RealTensorValue> &_fXf;
  MaterialProperty<RealTensorValue> &_gXg;
  MaterialProperty<RealTensorValue> &_sXs;
  MaterialProperty<RealTensorValue> &_fXs;
  MaterialProperty<RealTensorValue> &_nXn;
  MaterialProperty<RealTensorValue> &_rotationTensor;
};

#endif // FIBERSMATERIAL_H
