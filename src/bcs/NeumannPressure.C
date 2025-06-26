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

#include "NeumannPressure.h"
#include "Function.h"
#include "MooseMesh.h"

registerMooseObject("MECHApp", NeumannPressure);
template <>
InputParameters
validParams<NeumannPressure>()
{
    InputParameters params = validParams<IntegratedBC>();
    
    // Here we are adding a parameter that will be extracted from the input file
    // by the Parser
    params.addRequiredParam<Real>("pressure", "pressure");
    params.addRequiredParam<unsigned>("component", "component");
    params.addParam<FunctionName>("function",
                                  "The function that describes the pressure");
    
    params.addRequiredCoupledVar(
                                 "disp_x",
                                 "The gradient of this variable will be used as the velocity vector.");
    params.addCoupledVar(
                         "disp_y",
                         "The gradient of this variable will be used as the velocity vector.");
    params.addCoupledVar(
                         "disp_z",
                         "The gradient of this variable will be used as the velocity vector.");
    
    params.set<bool>("use_displaced_mesh") = false;
    
    return params;
}

NeumannPressure::NeumannPressure(const InputParameters &parameters)
: IntegratedBC(parameters), _pressure(getParam<Real>("pressure")),
_component(getParam<unsigned>("component")),
_function2(getFunction("function")),
_grad_disp_x(coupledGradient("disp_x")),
_grad_disp_y(_mesh.dimension() >= 2 ? coupledGradient("disp_y") : _grad_zero),
_grad_disp_z(_mesh.dimension() >= 3 ? coupledGradient("disp_z") : _grad_zero)

{
}

Real
NeumannPressure::computeQpResidual()
{
    
    RealTensorValue F(_grad_disp_x[_qp], _grad_disp_y[_qp], _grad_disp_z[_qp]);
    RealVectorValue normal;
    
    for (int i = 0; i < 3; ++i)
    {
        F(i, i) += 1.0;
        normal(i) = _normals[_qp](i);
    }
    
    Real J = F.det();
    
    RealTensorValue invFtr = F.inverse().transpose();
    
    RealVectorValue temp =  J * invFtr * normal;
    
    Real factor = _pressure * _function2.value(_t, _q_point[_qp]);
    
    return -factor * temp(_component) * _test[_i][_qp];
}

Real
NeumannPressure::computeQpJacobian()
{
    RealVectorValue zeros(0.0, 0.0, 0.0);
    
    RealTensorValue H(zeros, zeros, zeros);
    
    for (unsigned k = 0; k < 3; ++k)
        H(_component, k) = _grad_phi[_j][_qp](k);
    
    RealTensorValue F(_grad_disp_x[_qp], _grad_disp_y[_qp], _grad_disp_z[_qp]);
    RealVectorValue normal;
    
    for (int i = 0; i < 3; ++i)
    {
        F(i, i) += 1.0;
        normal(i) = _normals[_qp](i);
    }
    
    Real J = F.det();
    RealTensorValue invFtr = F.inverse().transpose();
    
    RealTensorValue TEMP =
    invFtr.contract(H) * invFtr - invFtr * H.transpose() * invFtr;
    RealVectorValue temp = J * TEMP * normal;
    
    Real factor = _pressure * _function2.value(_t, _q_point[_qp]);
    
    return -factor * temp(_component) * _test[_i][_qp];
}

Real
NeumannPressure::computeQpOffDiagJacobian(unsigned int jvar)
{
    if (jvar == 3)
        return 0.0;
    
    RealVectorValue zeros(0.0, 0.0, 0.0);
    
    RealTensorValue H(zeros, zeros, zeros);
    
    for (unsigned k = 0; k < 3; ++k)
        H(jvar, k) = _grad_phi[_j][_qp](k);
    
    RealTensorValue F(_grad_disp_x[_qp], _grad_disp_y[_qp], _grad_disp_z[_qp]);
    RealVectorValue normal;
    
    for (int i = 0; i < 3; ++i)
    {
        F(i, i) += 1.0;
        normal(i) = _normals[_qp](i);
    }
    
    Real J = F.det();
    RealTensorValue invFtr = F.inverse().transpose();
    
    RealTensorValue TEMP =
    (invFtr.contract(H)) * invFtr - (invFtr * (H.transpose())) * invFtr;
    RealVectorValue temp = J * TEMP * normal;
    
    Real factor = _pressure * _function2.value(_t, _q_point[_qp]);
    
    return -factor * temp(_component) * _test[_i][_qp];
}
