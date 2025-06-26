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
/*                Prepared by Marco Favino,                     */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/****************************************************************/


#include "ActiveStressFunction.h"

#include "Material.h"
#include "Function.h"

registerMooseObject("MECHApp", ActiveStressFunction);

template <>
InputParameters
validParams<ActiveStressFunction>()
{

  // inherit the parameters of the Kernels:
  InputParameters params = validParams<Kernel>();

  // specify for which component to solve:
  params.addRequiredParam<unsigned>("component", "component");

  // the driving function
  params.addParam<FunctionName>("function",
                                "The function that describes the pressure");
  return params;
}

ActiveStressFunction::ActiveStressFunction(const InputParameters &parameters)
    :

      // inherit the parameters of the Kernels:
      Kernel(parameters),

      // specify for which component to solve:
      _component(getParam<unsigned>("component")),

      // inherit some material properties:
      _fXf(getMaterialProperty<RealTensorValue>("fiberOuterFiber")),
      _F(getMaterialProperty<RealTensorValue>("deformationGradient")),

      // the driving function
      _activeFunction(getFunction("function"))
{
}

Real
ActiveStressFunction::computeQpResidual()
{

  // generate the test function
  V = RealTensorValue(zeros, zeros, zeros);

  for (unsigned k = 0; k < 3; ++k)
    V(_component, k) = _grad_test[_i][_qp](k);

  RealTensorValue Pa;

  // evaluation of the function at the current time for the quadrature point
  Real factor = _activeFunction.value(_t, _q_point[_qp]);

  RealTensorValue F = _F[_qp];

  // computation of the active tension
  Pa = factor * F * _fXf[_qp];

  // contraction of the active tension with the test function
  return Pa.contract(V);
}

Real
ActiveStressFunction::computeQpJacobian()
{

  // Computation of the diagonal entries of the Jacobian matrix

  H = RealTensorValue(zeros, zeros, zeros);
  V = RealTensorValue(zeros, zeros, zeros);

  for (unsigned k = 0; k < 3; ++k)
  {
    H(_component, k) = _grad_phi[_j][_qp](k);
  }

  for (unsigned k = 0; k < 3; ++k)
    V(_component, k) = _grad_test[_i][_qp](k);

  // evaluation of the function at the current time for the quadrature point
  Real factor = _activeFunction.value(_t, _q_point[_qp]);

  // computation of the linearization of the active tension
  RealTensorValue Pa = factor * H * _fXf[_qp];

  // contraction of the linearization of the active tension with the test
  // function
  return Pa.contract(V);
}

Real
ActiveStressFunction::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Computation of the off diagonal entries of the Jacobian matrix 
  // jvar is the variable w.r.t. which we are derivating (the considered coloumn). 

  if (jvar > 2)
    return 0.0;

  H = RealTensorValue(zeros, zeros, zeros);
  V = RealTensorValue(zeros, zeros, zeros);

  for (unsigned k = 0; k < 3; ++k)
  {
    H(jvar, k) = _grad_phi[_j][_qp](k);
  }

  for (unsigned k = 0; k < 3; ++k)
    V(_component, k) = _grad_test[_i][_qp](k);

  Real factor = _activeFunction.value(_t, _q_point[_qp]);

  // linearization of the active tension
  RealTensorValue Pa = factor * H * _fXf[_qp];

  // contraction of the linearized active tension with the test function
  return Pa.contract(V);
}
