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
/* Kernel to ensure conservation of momentum for incompressible */
/* hyperelastic materials. To be used only coupled with         */
/* MassBalance kernel.                                          */
/****************************************************************/

#include "IncompressibleStressDivergenceLibmeshTensor.h"
#include "Material.h"
#include "Assembly.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "SystemBase.h"

#include "libmesh/quadrature.h"
registerMooseObject("MECHApp", IncompressibleStressDivergenceLibmeshTensor);

template <>
InputParameters
validParams<IncompressibleStressDivergenceLibmeshTensor>()
{

  // inherit the parameters of the Kernels:
  InputParameters params = validParams<Kernel>();

  // specify for which component to solve:
  params.addRequiredParam<unsigned>("component", "component");

  // solve for the pressure in the massBalance Kernel implies pressure as coupled variable:
  params.addRequiredCoupledVar("pressure", "pressure");
  return params;
}

IncompressibleStressDivergenceLibmeshTensor::
    IncompressibleStressDivergenceLibmeshTensor(
        InputParameters const &parameters):

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

      // the pressure variable in coupled and non coupled copy:
      _p_var(coupled("pressure")), _p(coupledValue("pressure")),

      // inherit material properties:
      _J(getMaterialProperty<Real>("deformationDeterminant")),
      _invFtr(getMaterialProperty<RealTensorValue>("invFtr")),
      _P(getMaterialProperty<RealTensorValue>("firstPiolaKirchhofTensor")),

      // material pointer:
      _materialPointer(
          getMaterialProperty<ElasticMaterial *>("materialPointer"))
{

  // we suppose that the components of displacements are the first three:
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
IncompressibleStressDivergenceLibmeshTensor::computeQpResidual()
{

//if (_component==0 &&_qp==0 && _current_elem->id()==1)std::cout<<" _qp_elastic_ink = "<< _q_point[_qp] <<"  variable_elastic_ink = "<<_u[_qp]<<std::endl;
 

  // additive decomposition of Eq.
    RealTensorValue P = _P[_qp] + _p[_qp] * _J[_qp] * _invFtr[_qp];

  // generates the incompressible stress kernel of Eq.         for x_i=
  // _component
  // here the contraction with the gradient of the test function
  return P(_component, 0) * _grad_test[_i][_qp](0) +
         P(_component, 1) * _grad_test[_i][_qp](1) +
         P(_component, 2) * _grad_test[_i][_qp](2);
}

Real
IncompressibleStressDivergenceLibmeshTensor::computeQpJacobian()
{
  // Computation of the diagonal entries of the Jacobian matrix
    
//  if (_component==0 &&_qp==0 && _current_elem->id()==1)std::cout<<" _qp_elastic_ink_jac = "<< _q_point[_qp] <<"  variable_elastic_ink_jac = "<<_u[_qp]<<std::endl;

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

        // we need to add the contribution of the linearization of the term to
        // enforce incompressibility
        _stress_lin[j][qp] +=
            _p[qp] * _J[qp] * (_invFtr[qp].contract(H) * _invFtr[qp] -
                               _invFtr[qp] * H.transpose() * _invFtr[qp]);
      }
  }

  // here we contract the linearization of the Piola-Kirchhoff Stress Divergence
  // tensor with the test function.
  return _stress_lin[_j][_qp](_component, 0) * _grad_test[_i][_qp](0) +
         _stress_lin[_j][_qp](_component, 1) * _grad_test[_i][_qp](1) +
         _stress_lin[_j][_qp](_component, 2) * _grad_test[_i][_qp](2);
}

Real
IncompressibleStressDivergenceLibmeshTensor::computeQpOffDiagJacobian(
    unsigned int jvar)
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

          // we need to add the contribution of the linearization of the term to
          // enforce incompressibility
          _stress_lin[j][qp] +=
              _p[qp] * _J[qp] * (_invFtr[qp].contract(H) * _invFtr[qp] -
                                 _invFtr[qp] * H.transpose() * _invFtr[qp]);
        }
    }

    // here we contract the linearization of the Piola-Kirchhoff Stress
    // Divergence tensor with the gradient of the test function.
    return _stress_lin[_j][_qp](_component, 0) * _grad_test[_i][_qp](0) +
           _stress_lin[_j][_qp](_component, 1) * _grad_test[_i][_qp](1) +
           _stress_lin[_j][_qp](_component, 2) * _grad_test[_i][_qp](2);
  }

  // else for the variable pressure
  else
  {
    if (jvar == _p_var)
    {
      if (_i == 0 && _j == 0 && _qp == 0)
      {
        // initialize the tensor H
        H = RealTensorValue(zeros, zeros, zeros);

        // for all the components of the trial functions:
        for (int j = 0; j < _phi.size(); ++j)

          // for all the  quadrature points:
          for (int qp = 0; qp < _qrule->n_points(); ++qp)
          {
            // the derivative of disp_x, disp_y, disp_z w.r.t. p
            //(up to the trial function the evaluation is constant)
            _stress_lin[j][qp] = _phi[j][qp] * _J[qp] * _invFtr[qp];
          }
      }

      // here we contract the linearization of the Piola-Kirchhoff Stress
      // Divergence tensor with the test function.
      return _stress_lin[_j][_qp](_component, 0) * _grad_test[_i][_qp](0) +
             _stress_lin[_j][_qp](_component, 1) * _grad_test[_i][_qp](1) +
             _stress_lin[_j][_qp](_component, 2) * _grad_test[_i][_qp](2);
    }
    else
      return 0.0;
  }
}

void
IncompressibleStressDivergenceLibmeshTensor::computeJacobian()
{
  DenseMatrix<Number> &ke =
      _assembly.jacobianBlock(_var.number(), _var.number());
  _local_ke.resize(ke.m(), ke.n());
  _local_ke.zero();

  for (_i = 0; _i < _test.size(); _i++)
    for (_j = _i; _j < _phi.size(); _j++)
      for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        _local_ke(_i, _j) += _JxW[_qp] * _coord[_qp] * computeQpJacobian();

  for (_i = 0; _i < _test.size(); _i++)
    for (_j = 0; _j < _i; _j++)
      _local_ke(_i, _j) = _local_ke(_j, _i);

  ke += _local_ke;

  if (_has_diag_save_in)
  {
    unsigned int rows = ke.m();
    DenseVector<Number> diag(rows);
    for (unsigned int i = 0; i < rows; i++)
      diag(i) = _local_ke(i, i);

    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (unsigned int i = 0; i < _diag_save_in.size(); i++)
      _diag_save_in[i]->sys().solution().add_vector(
          diag, _diag_save_in[i]->dofIndices());
  }
}
