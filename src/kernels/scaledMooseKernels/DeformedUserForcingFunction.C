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

#include "DeformedUserForcingFunction.h"
#include "Function.h"

#include "MooseMesh.h"
#include "libmesh/quadrature.h"

registerMooseObject("MECHApp",DeformedUserForcingFunction);

template <>
InputParameters
validParams<DeformedUserForcingFunction>()
{

  // inherit the parameters of the Kernel:
  InputParameters params = validParams<Kernel>();

  
 // specify the function:
  params.addRequiredParam<FunctionName>("function", "The forcing function");
  
  return params;
}

DeformedUserForcingFunction::DeformedUserForcingFunction(
    const InputParameters &parameters):

     // inherit the parameters of the Kernel:
     Kernel(parameters),

      // We are doing elasticity from R^n to R^n,
      // hence the number of components is the same as the mesh dimension:
      // For example for _mesh.dimension=3 we have
      //_component = 0 for disp_x, _component = 1 for disp_y, _component = 2 for
      // disp_z
      _numComp(_mesh.dimension()),


      // specify the function:
      _func(getFunction("function")),

      // inherit some material properties. Problem solved the reference
      // configuration
      _J(getMaterialProperty<Real>("deformationDeterminant")),
      _invFtr(getMaterialProperty<RealTensorValue>("inverseTransposeStrain"))
{

  _J_lin = new Real *[64];
  for (int j = 0; j < 64; ++j)
    _J_lin[j] = new Real[64];
}

Real
DeformedUserForcingFunction::f()
{
  return _func.value(_t, _q_point[_qp]);
}

Real
DeformedUserForcingFunction::computeQpResidual()
{

  return -_J[_qp] *  _test[_i][_qp] * f();
}

Real
DeformedUserForcingFunction::computeQpJacobian()
{

  // computes the Jacobian in the reference configuration
  return 0.0;
}

Real
DeformedUserForcingFunction::computeQpOffDiagJacobian(unsigned int jvar)
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

          // evaluation of the linearization
          _J_lin[j][qp] = _J[qp] * _invFtr[qp].contract(H) * f();
        }
    }

    // here we contract the linearization of (J-1)
    // with the test function.
    return _J_lin[jvar][_qp] * _test[_i][_qp];
  }

  //   else

  return 0.0;
}
