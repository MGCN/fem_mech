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
/*           MECH - ICS Mechanical simulation framework         */
/*                Prepared by Maria Nestola                     */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/* Material to compute fibers direction.                        */
/****************************************************************/


#ifndef CompressibleHolzapfel_H
#define CompressibleHolzapfel_H

#include "ElasticMaterial.h"

// Forward declaration
class CompressibleHolzapfel;

template <>
InputParameters validParams<CompressibleHolzapfel>();


class CompressibleHolzapfel : public ElasticMaterial
{
public:
  CompressibleHolzapfel(const InputParameters &parameters);

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

  

  MaterialProperty<Real> &_c1;
  MaterialProperty<Real> &_c4;
  MaterialProperty<Real> &_c6;
  MaterialProperty<Real> &_c8;

 
  MaterialProperty<RealTensorValue> const &_fXf;
  MaterialProperty<RealTensorValue> const &_sXs;
  MaterialProperty<RealTensorValue> const &_fXs;
 
  MaterialProperty<RealTensorValue> &_C;
  MaterialProperty<RealTensorValue> &_invC;  


  MaterialProperty<RealTensorValue> &_S1;
  MaterialProperty<RealTensorValue> &_S4;
  MaterialProperty<RealTensorValue> &_S6;
  MaterialProperty<RealTensorValue> &_S8;

  MaterialProperty<RealTensorValue> &_S;

  
  RealTensorValue _S_lin;
  RealTensorValue _CLin;
  MaterialProperty<RealTensorValue> &_Pdev;



  
  void compute_C_lin(unsigned const &qp, RealTensorValue const &H);
};

#endif
