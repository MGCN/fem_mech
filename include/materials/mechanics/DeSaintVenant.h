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
/* Construction of Guccione Costa Material.                     */
/****************************************************************/

#ifndef DESAINTVENANT_H
#define DESAINTVENANT_H

#include "ElasticMaterial.h"

// Forward declaration
class DeSaintVenant;

template <>
InputParameters validParams<DeSaintVenant>();

/**
 * TensorMechanicsMaterial handles a fully anisotropic, single-crystal
 * material's elastic
 * constants.  It takes all 21 independent stiffness tensor inputs, or only 9,
 * depending on the
 * boolean input value given.  This can be extended or simplified to specify
 * HCP, monoclinic,
 * cubic, etc as needed.
 */
class DeSaintVenant : public ElasticMaterial
{
public:
  DeSaintVenant(const InputParameters &parameters);

  virtual void evaluate_stress_lin(unsigned const &qp, RealTensorValue const &H,
                                   RealTensorValue &stressLin);

  virtual void computeQpPropertiesDerived();

protected:
  //    virtual void computeQpProperties();

  Real mu;
  Real kappa;
  Real lambda;
  Real epsilon;

  MaterialProperty<RealTensorValue> &_E;

  MaterialProperty<RealTensorValue> &_C;

  MaterialProperty<RealTensorValue> &_S;
  
  MaterialProperty<RealTensorValue> &_Pdev;
};

#endif
