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

#ifndef ACTIVESTRESSFUNCTION_TIMEACTIVATED_H
#define ACTIVESTRESSFUNCTION_TIMEACTIVATED_H

#include "Kernel.h"
#include "ElasticMaterial.h"

class ActiveStressFunction_timeactivated;
class Function;

template <>
InputParameters validParams<ActiveStressFunction_timeactivated>();

class ActiveStressFunction_timeactivated : public Kernel
{
public:
  ActiveStressFunction_timeactivated(const InputParameters &parameters);

  unsigned int _component;

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const MaterialProperty<RealTensorValue> &_fXf;

  const MaterialProperty<RealTensorValue> &_F;

  RealVectorValue zeros;

  RealTensorValue H;
  RealTensorValue V;

  //const Function &_activeFunction_modulus;
 Real _activeFunction_modulus;
  const VariableValue &_time_activation_map;
  bool _time_unit;
};
#endif
