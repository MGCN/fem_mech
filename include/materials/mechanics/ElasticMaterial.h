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
/*                Prepared by Marco Favino,                     */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/* Construction of Guccione Costa Material.                     */
/****************************************************************/

#ifndef ELASTICMATERIAL_H
#define ELASTICMATERIAL_H

#include "Material.h"

// typedef RealTensorValue RTV[4];

class ElasticMaterial;

template <>
InputParameters validParams<ElasticMaterial>();

/**
 * TensorMechanicsMaterial handles a fully anisotropic, single-crystal
 * material's elastic
 * constants.  It takes all 21 independent stiffness tensor inputs, or only 9,
 * depending on the
 * boolean input value given.  This can be extended or simplified to specify
 * HCP, monoclinic,
 * cubic, etc as needed.
 */
class ElasticMaterial : public Material
{
public:
  ElasticMaterial(InputParameters const &parameters);

  void computeQpProperties();

  virtual void computeQpPropertiesDerived() = 0;

  virtual void evaluate_stress_lin(unsigned const &qp, RealTensorValue const &H,
                                   RealTensorValue &stressLin) = 0;

  virtual void initQpStatefulProperties();

protected:
    
    int const _dim;
    MaterialProperty<Real> &_J;
    MaterialProperty<RealTensorValue> &_Ce;
    MaterialProperty<RealTensorValue> &_Se;
    MaterialProperty<Real> &_pressure_material;
    
    
  const VariableGradient &_grad_disp_x;
  const VariableGradient &_grad_disp_y;
  const VariableGradient &_grad_disp_z;

  bool const _has_pressure;
  const VariableValue &_pressure;
    
    
    
  /// Material property base name to allow for multiple TensorMechanicsMaterial
  /// to coexist in the same simulation
  //    std::string _base_name;

  MaterialProperty<RealTensorValue> &_U;
  MaterialProperty<RealTensorValue> &_F;
  MaterialProperty<RealTensorValue> &_invFtr;
  MaterialProperty<RealTensorValue> &_P;
  MaterialProperty<RealTensorValue> &_stress;

  bool _store_variable_older;

  RealTensorValue _Id;

  //    MaterialProperty<RealTensorValue> & _applied_deformation_tensor;

  MaterialProperty<ElasticMaterial *> &materialPointer;

  //  bool _has_inelastic_deformation;
};

#endif
