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

#ifndef IncompressibleMooneyRivlinTP_H
#define IncompressibleMooneyRivlinTP_H

#include "ElasticMaterial.h"

// Forward declaration
class IncompressibleMooneyRivlinTP;

template <>
InputParameters validParams<IncompressibleMooneyRivlinTP>();

class IncompressibleMooneyRivlinTP : public ElasticMaterial
{
public:
  IncompressibleMooneyRivlinTP(const InputParameters &parameters);

  virtual void evaluate_stress_lin(unsigned const &qp, RealTensorValue const &H,
                                   RealTensorValue &stressLin);

                                   
  virtual void computeQpPropertiesDerived();

  void compute_E_prime(unsigned const &qp, RealTensorValue const &H);

protected:
  //    virtual void computeQpProperties();

  Real _c10;
  Real _c20;
  Real _c11;
  MaterialProperty<Real> &_I1_bar;
  MaterialProperty<Real> &_I2_bar;
  RealTensorValue Clin;
  RealTensorValue elin;
  RealTensorValue _Sbar_lin;
  RealTensorValue _Siso_lin;


  MaterialProperty<RealTensorValue> &_C;
  MaterialProperty<RealTensorValue> &_C_bar;
  MaterialProperty<RealTensorValue> &_invC;

  MaterialProperty<Real> &_J23;

  MaterialProperty<RealTensorValue> &_Sbar;

  MaterialProperty<RealTensorValue> &_Siso;


};

#endif