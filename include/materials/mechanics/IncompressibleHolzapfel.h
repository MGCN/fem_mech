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


#ifndef IncompressibleHolzapfel_H
#define IncompressibleHolzapfel_H

#include "ElasticMaterial.h"

// Forward declaration
class IncompressibleHolzapfel;

template <>
InputParameters validParams<IncompressibleHolzapfel>();


class IncompressibleHolzapfel : public ElasticMaterial
{
public:
  IncompressibleHolzapfel(const InputParameters &parameters);

  virtual void evaluate_stress_lin(unsigned const &qp, RealTensorValue const &H,
                                   RealTensorValue &stressLin);
  

  virtual void computeQpPropertiesDerived();

protected:
  //    virtual void computeQpProperties();

  Real mu;
  Real k_1f, k_2f;
  Real k_1s, k_2s;
  Real k_1fs, k_2fs;

  
  MaterialProperty<Real> &_I4;
  MaterialProperty<Real> &_I6;
  MaterialProperty<Real> &_I8;
  MaterialProperty<Real> &_I4_bar;
  MaterialProperty<Real> &_I6_bar;
  MaterialProperty<Real> &_I8_bar;
  

  MaterialProperty<Real> &_c1;
  MaterialProperty<Real> &_c4;
  MaterialProperty<Real> &_c6;
  MaterialProperty<Real> &_c8;

 
  MaterialProperty<RealTensorValue> const &_fXf;
  MaterialProperty<RealTensorValue> const &_gXg;
  MaterialProperty<RealTensorValue> const &_sXs;
  MaterialProperty<RealTensorValue> const &_fXs;
 
  MaterialProperty<RealTensorValue> &_C;
  MaterialProperty<RealTensorValue> &_invC;  

  MaterialProperty<Real> &_J23;
    
  MaterialProperty<RealTensorValue> &_S1bar;
  MaterialProperty<RealTensorValue> &_S4bar;
  MaterialProperty<RealTensorValue> &_S6bar;
  MaterialProperty<RealTensorValue> &_S8bar;

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
