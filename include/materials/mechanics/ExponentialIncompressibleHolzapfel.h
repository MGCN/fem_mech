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


#ifndef ExponentialIncompressibleHolzapfel_H
#define ExponentialIncompressibleHolzapfel_H

#include "ElasticMaterial.h"

// Forward declaration
class ExponentialIncompressibleHolzapfel;

template <>
InputParameters validParams<ExponentialIncompressibleHolzapfel>();


class ExponentialIncompressibleHolzapfel : public ElasticMaterial
{
public:
ExponentialIncompressibleHolzapfel(const InputParameters &parameters);

  virtual void evaluate_stress_lin(unsigned const &qp, RealTensorValue const &H,
                                   RealTensorValue &stressLin);
  

  virtual void computeQpPropertiesDerived();

protected:
  //    virtual void computeQpProperties();

  Real c01, c10;
  Real k_1f, k_2f;

  
  MaterialProperty<Real> &_I1;
  MaterialProperty<Real> &_I4;

  MaterialProperty<Real> &_I1_bar;
  MaterialProperty<Real> &_I4_bar;


  MaterialProperty<Real> &_c1;
  MaterialProperty<Real> &_c4_iso;
  MaterialProperty<Real> &_c4_aniso;
  MaterialProperty<Real> & kappa;


 
  MaterialProperty<RealTensorValue> const &_fXf;

 
  MaterialProperty<RealTensorValue> &_C;
  MaterialProperty<RealTensorValue> &_invC;  

  MaterialProperty<Real> &_J23;
    
  MaterialProperty<RealTensorValue> &_S1bar;
  MaterialProperty<RealTensorValue> &_S4bar;


  MaterialProperty<RealTensorValue> &_Sbar;

  MaterialProperty<RealTensorValue> &_Siso;
  
  RealTensorValue _Sbar_lin;
  RealTensorValue _Siso_lin;
  RealTensorValue  Cbar_lin; 
  RealTensorValue _CLin;
  RealTensorValue elin;
  MaterialProperty<RealTensorValue> &_Pdev;



  
  void compute_C_lin(unsigned const &qp, RealTensorValue const &H);
};

#endif
