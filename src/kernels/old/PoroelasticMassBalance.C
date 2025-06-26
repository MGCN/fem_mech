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
/* Kernel to ensure poroelastic mass balance.                   */
/****************************************************************/


#include "PoroelasticMassBalance.h"
#include "MooseMesh.h"

registerMooseObject("MECHApp", PoroelasticMassBalance);

template <>
InputParameters
validParams<PoroelasticMassBalance>()
{

  // inherit the parameters of the TimeKernel:
  InputParameters params = validParams<TimeKernel>();

  // specify for which component to solve:
  params.addRequiredParam<unsigned>("component", "component");
  return params;
}

PoroelasticMassBalance::PoroelasticMassBalance(
    InputParameters const &parameters)
    : 
   // inherit the parameters of the Kernels:
    TimeKernel(parameters),

  // specify for which component to solve:
     _component(getParam<unsigned>("component")),
    

      // We are doing elasticity from R^n to R^n,
      // hence the number of components is the same as the mesh dimension:
      // For example for _mesh.dimension=3 we have
      //_component = 0 for disp_x, _component = 1 for disp_y, _component = 2 for
      // disp_z
      _numComp(_mesh.dimension()),
      

    // inherit material properties:
      _J(getMaterialProperty<Real>("deformationDeterminant")),
      _J_old(getMaterialPropertyOld<Real>("deformationDeterminant")),
      _invFtr(getMaterialProperty<RealTensorValue>("invFtr"))
{

  // we suppose that the components of displacements are the first three.
  for (int i = 0; i < 3; ++i)
    zeros(i) = 0.0;
}

Real
PoroelasticMassBalance::computeQpResidual()
{
  return (_J[_qp] - _J_old[_qp]) * _test[_i][_qp];
}

Real
PoroelasticMassBalance::computeQpJacobian()
{
  return 0.0;
}

Real
PoroelasticMassBalance::computeQpOffDiagJacobian(unsigned int jvar)
{

  // Computation of the off diagonal entries of the Jacobian matrix 
  // jvar is the variable w.r.t. which we are derivating (the considered coloumn). 
  // for all the variables different from the pressure (the current one is jvar)
  // for example for mesh_dimension=3

  if (jvar < _numComp)
  {
    H = RealTensorValue(zeros, zeros, zeros);

    for (unsigned k = 0; k < _numComp; ++k)
    {
      H(jvar, k) = _grad_phi[_j][_qp](k);
    }

    return _J[_qp] * (_invFtr[_qp].contract(H)) * _test[_i][_qp];
  }
  else
    return 0.0;
}
