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
/*                Prepared by Maria Nestola,                    */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/****************************************************************/

#ifndef INCOMPRESSIBLENEOHOOKEAN_H
#define INCOMPRESSIBLENEOHOOKEAN_H

#include "ElasticMaterial.h"

// Forward declaration
class IncompressibleNeoHookean;

template <>
InputParameters validParams<IncompressibleNeoHookean>();

/**
 * TensorMechanicsMaterial handles a fully anisotropic, single-crystal
 * material's elastic
 * constants.  It takes all 21 independent stiffness tensor inputs, or only 9,
 * depending on the
 * boolean input value given.  This can be extended or simplified to specify
 * HCP, monoclinic,
 * cubic, etc as needed.
 */
class IncompressibleNeoHookean : public ElasticMaterial
{
public:
  IncompressibleNeoHookean(const InputParameters &parameters);

  virtual void evaluate_stress_lin(unsigned const &qp, RealTensorValue const &H,
                                   RealTensorValue &stressLin);

  virtual void computeQpPropertiesDerived();

  void compute_E_prime(unsigned const &qp, RealTensorValue const &H);

protected:
  //    virtual void computeQpProperties();

  Real mu;

  RealTensorValue Clin;
  RealTensorValue _Siso_lin;

  MaterialProperty<RealTensorValue> &_C;
  MaterialProperty<RealTensorValue> &_invC;

  MaterialProperty<Real> &_J23;

  MaterialProperty<RealTensorValue> &_Sbar;

  MaterialProperty<RealTensorValue> &_Siso;
    
  MaterialProperty<RealTensorValue> &_Pdev;
  
 /// const Function &_function2;    
};

#endif
