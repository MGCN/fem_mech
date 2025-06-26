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
/* Kernel to ensure conservation of momentum for compressible   */
/* hyperelastic materials.                                      */
/****************************************************************/

#include "StressDivergenceLibmeshTensor.h"
#include "Material.h"
#include "MooseMesh.h"
#include "libmesh/quadrature.h"

registerMooseObject("MECHApp", StressDivergenceLibmeshTensor);

template <>
InputParameters
validParams<StressDivergenceLibmeshTensor>()
{

  // inherit the parameters of the Kernels:
  InputParameters params = validParams<Kernel>();

 // specify for which component to solve:
  params.addRequiredParam<unsigned>("component", "component");

  return params;
}

StressDivergenceLibmeshTensor::StressDivergenceLibmeshTensor(
    InputParameters const &parameters)
    : // inherit the parameters of the Kernels:
      Kernel(parameters),

      // specify for which component to solve:
      _component(getParam<unsigned>("component")),

      // inherit some material properties:
      _P(getMaterialProperty<RealTensorValue>("firstPiolaKirchhofTensor")),

      // material pointer:
      _materialPointer(
          getMaterialProperty<ElasticMaterial *>("materialPointer")),

      // We are doing elasticity from R^n to R^n,
      // hence the number of components is the same as the mesh dimension:
      // For example for _mesh.dimension=3 we have
      //_component = 0 for disp_x, _component = 1 for disp_y, _component = 2 for
      // disp_z
      _numComp(_mesh.dimension())
{
  // we suppose that the components of displacements are the first three.
  if (_component >= _numComp)
  {
    mooseError("Error in IncompressibleStressDivergenceLibmeshTensor: one of "
               "the component is larger than mesh dimension.");
  }

  for (int i = 0; i < 3; ++i)
    zeros(i) = 0.0;
  _stress_lin = new RealTensorValue *[64] /*[_phi.size()]*/;

  for (int j = 0; j < 64 /*_phi.size()*/; ++j)
    _stress_lin[j] = new RealTensorValue[64] /*[_qrule->n_points()]*/;


}

Real
StressDivergenceLibmeshTensor::computeQpResidual()
{

  // generates the compressible stress kernel of Eq. for x_i= _component
  // here is performed the contraction with the gradient if the test functin
  return _P[_qp](_component, 0) * _grad_test[_i][_qp](0) +
         _P[_qp](_component, 1) * _grad_test[_i][_qp](1) +
         _P[_qp](_component, 2) * _grad_test[_i][_qp](2);
}

Real
StressDivergenceLibmeshTensor::computeQpJacobian()
{

  // Computation of the diagonal entries of the Jacobian matrix
  if (_i == 0 && _j == 0 && _qp == 0)
  {
    // initialize the trial tensor H
    H = RealTensorValue(zeros, zeros, zeros);

    // for all the components of the trial functions
    for (int j = 0; j < _phi.size(); ++j)

      // for all the quadrature points
      for (int qp = 0; qp < _qrule->n_points(); ++qp)
      {

        // in H we store the gradient of trial function: in this way j becomes
        // _component for H
        for (unsigned k = 0; k < _numComp; ++k)
        {
          H(_component, k) = _grad_phi[j][qp](k);
        }

        // evaluation of the linearization of the stress given by the material
        (*(_materialPointer[_qp]))
            .evaluate_stress_lin(qp, H, _stress_lin[j][qp]);
      }
  }

  // here we contract the linearization of the Piola-Kirchhoff Stress Divergence
  // tensor with the test function.
  return _stress_lin[_j][_qp](_component, 0) * _grad_test[_i][_qp](0) +
         _stress_lin[_j][_qp](_component, 1) * _grad_test[_i][_qp](1) +
         _stress_lin[_j][_qp](_component, 2) * _grad_test[_i][_qp](2);
}

Real
StressDivergenceLibmeshTensor::computeQpOffDiagJacobian(unsigned int jvar)
{ 
  // Computation of the off diagonal entries of the Jacobian matrix 
  // jvar is the variable w.r.t. which we are derivating (the considered coloumn). 
  // for all the variables different from the pressure (the current one is jvar)
  // for example for mesh_dimension=3


// if (jvar==3)
//   return 0.0;
if (jvar < _numComp){
  if (_i == 0 && _j == 0 && _qp == 0)
  {

    // for all the components of the trial functions:
    H = RealTensorValue(zeros, zeros, zeros);

    // for all the components of the trial functions:
    for (int j = 0; j < _phi.size(); ++j)

      // for all the  quadrature points:
      for (int qp = 0; qp < _qrule->n_points(); ++qp)
      {

        // in H we store the gradient of trial function: in this way j becomes
        // _jvar for H
        for (unsigned k = 0; k < _numComp; ++k)
        {
          H(jvar, k) = _grad_phi[j][qp](k);
        }

        // evaluation of the linearization of the stress given by the material
        (*(_materialPointer[_qp]))
            .evaluate_stress_lin(qp, H, _stress_lin[j][qp]);
      }
  }
}else
{ 
  return 0.0; 
}

  // here we contract the linearization of the Piola-Kirchhoff Stress
  // Divergence tensor with the gradient of the test function.
  return _stress_lin[_j][_qp](_component, 0) * _grad_test[_i][_qp](0) +
         _stress_lin[_j][_qp](_component, 1) * _grad_test[_i][_qp](1) +
         _stress_lin[_j][_qp](_component, 2) * _grad_test[_i][_qp](2);
}
