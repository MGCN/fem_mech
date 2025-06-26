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

#ifndef COMPRESSIBLEGUCCIONECOSTA_H
#define COMPRESSIBLEGUCCIONECOSTA_H

#include "ElasticMaterial.h"

// Forward declaration
class CompressibleGuccioneCosta;

template <>
InputParameters validParams<CompressibleGuccioneCosta>();


class CompressibleGuccioneCosta : public ElasticMaterial
{
public:
  CompressibleGuccioneCosta(InputParameters const &parameters);

  virtual void computeQpPropertiesDerived();

  virtual void evaluate_stress_lin(unsigned const &qp, RealTensorValue const &H,
                                   RealTensorValue &stressLin);

  void compute_E_prime(unsigned const &qp, RealTensorValue const &H);

protected:
  //    virtual void computeQpProperties();

  Real mu;
  Real bf;
  Real bfs;
  Real bt;
  RealTensorValue _Q;

  RealTensorValue _EH;

  MaterialProperty<RealTensorValue> &_E;

  MaterialProperty<RealTensorValue> &_C;

  MaterialProperty<RealTensorValue> &_EtimesQ;

  MaterialProperty<RealTensorValue> &_S;

  MaterialProperty<Real> &_stiffening;
};

#endif