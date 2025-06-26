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

#ifndef GrowthMaterial_H
#define GrowthMaterial_H

#include "ElasticMaterial.h"

// Forward declaration
class GrowthMaterial;

template <>
InputParameters validParams<GrowthMaterial>();

class GrowthMaterial : public ElasticMaterial
{
public:
 GrowthMaterial(const InputParameters &parameters);

  virtual void evaluate_stress_lin(unsigned const &qp, RealTensorValue const &H,
                                   RealTensorValue &stressLin);

                                   
  virtual void computeQpPropertiesDerived();

  //void compute_E_prime(unsigned const &qp, RealTensorValue const &H);

protected:
  //    virtual void computeQpProperties();

  //Real _c10;
  Real _mu;
  Real _lambda;

  const VariableValue & _theta;
  
  // RealTensorValue Clin;

  MaterialProperty<RealTensorValue> &_C;

  MaterialProperty<RealTensorValue> &_invC;

  MaterialProperty<RealTensorValue> &_S;

  // MaterialProperty<RealTensorValue> &_Ce;

  MaterialProperty<RealTensorValue> &_invCe;

  // MaterialProperty<RealTensorValue> &_Se;

  MaterialProperty<RealTensorValue> &_Fg;

  MaterialProperty<RealTensorValue> &_Fgt;  

  MaterialProperty<RealTensorValue> &_invFg;

  MaterialProperty<RealTensorValue> &_invFgt;

};

#endif