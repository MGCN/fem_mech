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
/* Kernel to compute the time derivative of come variable.      */
/* The solution is computed on the deformed domain.             */
/****************************************************************/
#include "DeformedTimeDerivative.h"
#include "Material.h"

#include "MooseMesh.h"
#include "libmesh/quadrature.h"

registerMooseObject("MECHApp", DeformedTimeDerivative);
template <>
InputParameters
validParams<DeformedTimeDerivative>()
{

  // inherit the parameters of the TimeDerivative:
  InputParameters params = validParams<TimeDerivative>();

   // inherit the parameters of the Kernel:
  params += validParams<Kernel>();

  params.addRequiredParam<Real>(
      "surface_to_volume",
      "Surface-to-volume ratio in 1/cm. Suggestion: 1400.0 ");
  // specify the capacitance:
  params.addRequiredParam<Real>(
      "capacitance",
      "Capacitance of membrane in microF/mm^2.");
  return params;
}

DeformedTimeDerivative::DeformedTimeDerivative(
    const InputParameters &parameters):

      // inherit the parameters of the TimeDerivative:
      TimeDerivative(parameters),

      // We are doing elasticity from R^n to R^n,
      // hence the number of components is the same as the mesh dimension:
      // For example for _mesh.dimension=3 we have
      //_component = 0 for disp_x, _component = 1 for disp_y, _component = 2 for
      // disp_z
      _numComp(_mesh.dimension()),

      // inherit some material properties:
      _J(getMaterialProperty<Real>("deformationDeterminant")),
      _invFtr(getMaterialProperty<RealTensorValue>("inverseTransposeStrain")),
      _surface_to_volume(getParam<Real>("surface_to_volume")),
      // specify the capacitance:
      _capacitance(getParam<Real>("capacitance"))
{

  // we suppose that the components of displacements are the first three.
  for (int i = 0; i < 3; ++i)
    zeros(i) = 0.0;

  _J_lin = new Real *[64];
  for (int j = 0; j < 64; ++j)
    _J_lin[j] = new Real[64];
}

Real
DeformedTimeDerivative::computeQpResidual()
{

  // computes the residual in the reference configuration through
  // TimeDerivative::computeQpResidual()
  return _J[_qp] *_surface_to_volume * _capacitance * TimeDerivative::computeQpResidual();
}

Real
DeformedTimeDerivative::computeQpJacobian()
{

  // computes the jacobian in the reference configuration through
  // TimeDerivative::computeQpJacobian()
  return _J[_qp] * _surface_to_volume *_capacitance * TimeDerivative::computeQpJacobian();
}

Real
DeformedTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
    // Computation of the off diagonal entries of the Jacobian matrix 
  // jvar is the variable w.r.t. which we are derivating (the considered coloumn). 
  // for all the variables different from the pressure (the current one is jvar)
  // for example for mesh_dimension=3

  if (jvar < _numComp)
  {
    if (_i == 0 && _j == 0 && _qp == 0)
    {
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
          _J_lin[j][qp] = _J[_qp] * _invFtr[_qp].contract(H) * _surface_to_volume * _capacitance *
                          TimeDerivative::computeQpResidual();
        }
    }

    // here we contract the linearization of (J-1)
    // with the test function.
    return _J_lin[_j][_qp] * _test[_i][_qp];
  }

  else

    return 0.0;
}
