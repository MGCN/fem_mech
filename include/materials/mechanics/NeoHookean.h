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
/*                Prepared by Maria Nestola,                    */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/****************************************************************/



#ifndef NEOHOOKEAN_H
#define NEOHOOKEAN_H

#include "ElasticMaterial.h"
#include "MooseMesh.h"

// Forward declaration
class NeoHookean;

template <>
InputParameters validParams<NeoHookean>();

/**
 * TensorMechanicsMaterial handles a fully anisotropic, single-crystal
 * material's elastic
 * constants.  It takes all 21 independent stiffness tensor inputs, or only 9,
 * depending on the
 * boolean ß value given.  This can be extended or simplified to specify HCP,
 * monoclinic,
 * cubic, etc as needed.
 */
class NeoHookean : public ElasticMaterial
{
public:
  NeoHookean(const InputParameters &parameters);

  virtual void evaluate_stress_lin(unsigned const &qp, RealTensorValue const &H,
                                   RealTensorValue &stressLin);

  virtual void computeQpPropertiesDerived();

protected:
  Real mu;
  Real lambda;
  MaterialProperty<RealTensorValue> &_C;
  MaterialProperty<RealTensorValue> &_invC;
};

#endif
