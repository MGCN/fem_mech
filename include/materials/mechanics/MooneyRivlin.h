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

#ifndef MooneyRivlin_H
#define MooneyRivlin_H

#include "ElasticMaterial.h"

// Forward declaration
class MooneyRivlin;

template <>
InputParameters validParams<MooneyRivlin>();

class MooneyRivlin : public ElasticMaterial
{
public:
 MooneyRivlin(const InputParameters &parameters);

  virtual void evaluate_stress_lin(unsigned const &qp, RealTensorValue const &H,
                                   RealTensorValue &stressLin);

                                   
  virtual void computeQpPropertiesDerived();

  void compute_E_prime(unsigned const &qp, RealTensorValue const &H);

protected:
  //    virtual void computeQpProperties();

  Real _c10;
  Real _c20;
  Real _lambda;
  
  RealTensorValue Clin;

  MaterialProperty<RealTensorValue> &_C;
  MaterialProperty<RealTensorValue> &_invC;

};

#endif