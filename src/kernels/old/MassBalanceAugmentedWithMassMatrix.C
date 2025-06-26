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
/* Kernel to ensure the mass balance. To be used only coupled   */
/* with IncompressibleStressDivergenceLibmeshTensor             */
/****************************************************************/

#include "MassBalanceAugmentedWithMassMatrix.h"
#include "MooseMesh.h"

#include "libmesh/quadrature.h"
registerMooseObject("MECHApp", MassBalanceAugmentedWithMassMatrix);
template <>
InputParameters
validParams<MassBalanceAugmentedWithMassMatrix>()
{

  // inherit the parameters of the Kernels:
  InputParameters params = validParams<Kernel>();

  // specify for which component to solve:
  params.addRequiredParam<unsigned>("component", "component");

  return params;
}

MassBalanceAugmentedWithMassMatrix::MassBalanceAugmentedWithMassMatrix(
    InputParameters const & parameters)
  :

    // inherit the parameters of the Kernels:
    Kernel(parameters),

    // we need to specify for which component we want to solve:
    _component(getParam<unsigned>("component")),

    // We are doing elasticity from R^n to R^n,
    // hence the number of components is the same as the mesh dimension:
    // For example for _mesh.dimension=3 we have
    //_component = 0 for disp_x, _component = 1 for disp_y, _component = 2 for
    // disp_z

    _numComp(_mesh.dimension()),

    // inherit some material properties:
    _J(getMaterialProperty<Real>("deformationDeterminant")),
    _invFtr(getMaterialProperty<RealTensorValue>("invFtr"))
{

  // we suppose that the components of displacements are the first three.
  for (int i = 0; i < 3; ++i)
    zeros(i) = 0.0;

  _J_lin = new Real *[64] /*[_phi.size()]*/;
  for (int j = 0; j < 64 /*_phi.size()*/; ++j)
    _J_lin[j] = new Real[64] /*[_qrule->n_points()]*/;
}

Real
MassBalanceAugmentedWithMassMatrix::computeQpResidual()
{
  return (_J[_qp] - 1.0) * _test[_i][_qp];
}

Real
MassBalanceAugmentedWithMassMatrix::computeQpJacobian()
{
  // the derivative of (J-1) w.r.t p is 0 but here we return the mass matrix
  return _test[_i][_qp] * _phi[_j][_qp];
}

Real
MassBalanceAugmentedWithMassMatrix::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Computation of the off diagonal entries of the Jacobian matrix
  // jvar is the variable w.r.t. which we are derivating (the considered
  // coloumn).
  // In this case in three-dimension, e.g. the derivative of (J-1) w.r.t disp_x,
  // disp_y, disp_z

  if (jvar < _numComp) {
    if (_i == 0 && _j == 0 && _qp == 0) {
      H = RealTensorValue(zeros, zeros, zeros);

      // for all the components of the trial functions:
      for (int j = 0; j < _phi.size(); ++j)

        // for all the quadrature points:
        for (int qp = 0; qp < _qrule->n_points(); ++qp)
        {

          // in H we store the gradient of trial function: in this way j becomes
          // _jvar for H
          for (unsigned k = 0; k < _numComp; ++k)
          {
            H(jvar, k) = _grad_phi[j][qp](k);
          }

          // evaluation of the linearization of (J-1)
          _J_lin[j][qp] = _J[qp] * _invFtr[qp].contract(H);
        }
    }

    // here we contract the linearization of (J-1)
    // with the test function.
    return _J_lin[_j][_qp] * _test[_i][_qp];
  }
  else
    return 0.0;
}
