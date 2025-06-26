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

#ifndef GUCCIONECOSTAROT_H
#define GUCCIONECOSTAROT_H

#include "ElasticMaterial.h"

// Forward declaration
class GuccioneCostaRot;

template <>
InputParameters validParams<GuccioneCostaRot>();

/**
 * TensorMechanicsMaterial handles a fully anisotropic, single-crystal
 * material's elastic
 * constants.  It takes all 21 independent stiffness tensor inputs, or only 9,
 * depending on the
 * boolean input value given.  This can be extended or simplified to specify
 * HCP, monoclinic,
 * cubic, etc as needed.
 */
class GuccioneCostaRot : public ElasticMaterial
{
public:
  GuccioneCostaRot(const InputParameters &parameters);

  virtual void computeQpPropertiesDerived();

  virtual void evaluate_stress_lin(unsigned const &qp, RealTensorValue const &H,
                                   RealTensorValue &stressLin);

  void compute_Ebar_prime(unsigned const &qp, RealTensorValue const &H);

protected:
  //    virtual void computeQpProperties();

  Real mu;
  Real bf;
  Real bfs;
  Real bt;
  RealTensorValue _Q;
  RealTensorValue _CH;

  MaterialProperty<RealTensorValue> &_Ebar;

  MaterialProperty<RealTensorValue> &_EbarRot;

  MaterialProperty<RealTensorValue> &_C;

  MaterialProperty<RealTensorValue> &_invC;

  MaterialProperty<RealTensorValue> &_EtimesQ;

  MaterialProperty<RealTensorValue> &_SbarRot;

  MaterialProperty<RealTensorValue> &_Sbar;

  MaterialProperty<RealTensorValue> &_S;
  MaterialProperty<RealTensorValue> &_Pdev;

  MaterialProperty<RealTensorValue> const &_R;

  MaterialProperty<Real> &_J23;

  MaterialProperty<Real> &_stiffening;

  RealTensorValue _EbarH;

  RealTensorValue _EbarH_Rot;

  RealTensorValue _EH_Times_Q;

  RealTensorValue _SbarH;

  //  void computeEPrime(RealTensorValue const &input, RealTensorValue &output);
  //
  //  void computeESecond(RealTensorValue const &input1,
  //                      RealTensorValue const &input2, RealTensorValue
  //                      &output);
};

#endif
