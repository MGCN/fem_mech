/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*        MECH - ICS Imm Bound  simulation framework            */
/*                Prepared by Maria Nestola                     */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/*         Kernel for the Nearly Incompressibility              */
/****************************************************************/
#include "RobinCondition.h"
#include "MooseMesh.h"
#include "Material.h"
#include "ElasticMaterial.h"
#include "libmesh/quadrature.h"
#include "Function.h"

registerMooseObject("MECHApp", RobinCondition);
template<>
InputParameters validParams<RobinCondition>()
{
    InputParameters params = validParams<IntegratedBC>();
   
    params.addParam<Real>("alpha_disp", 1.0,
                          "alpha_disp");
    params.addParam<unsigned int>("component", "component");
    params.addRequiredParam<FunctionName>("function", "The function.");
  
    return params;
}

RobinCondition::RobinCondition(const InputParameters & parameters) :
IntegratedBC(parameters),
_component(getParam<unsigned int>("component")),
_alpha_disp(getParam<Real>("alpha_disp")),
// _alpha_stress(getParam<Real>("alpha_stress")),
_pext(getFunction("function"))


{
}

Real
RobinCondition::computeQpResidual()
{
    
  return _alpha_disp * _u[_qp] * _test[_i][_qp] - _pext.value(_t, _q_point[_qp]) * _normals[_qp](_component) * _test[_i][_qp] ;
}

Real
RobinCondition::computeQpJacobian()
{
  return _alpha_disp * _phi[_j][_qp] * _test[_i][_qp];
}





