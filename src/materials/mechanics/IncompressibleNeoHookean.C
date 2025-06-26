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
/*                Prepared by Maria Nestola,                    */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/****************************************************************/

#include "IncompressibleNeoHookean.h"
#include "MooseMesh.h"
#include "Function.h"
registerMooseObject("MECHApp", IncompressibleNeoHookean);

template <>
InputParameters
validParams<IncompressibleNeoHookean>()
{
  InputParameters params = validParams<ElasticMaterial>();
  params.addRequiredParam<Real>("mu", "mu");
 // params.addParam<FunctionName>("function",
  //                                "The function that describes the pressure");
    
  return params;
}
 
IncompressibleNeoHookean::IncompressibleNeoHookean(const InputParameters &parameters) :
    ElasticMaterial(parameters),
    mu(getParam<Real>("mu")),
    _C(declareProperty<RealTensorValue>("CauchyGreen")),
    _invC(declareProperty<RealTensorValue>("inverseCauchyGreen")),
    _J23(declareProperty<Real>("J23")),
    _Sbar(declareProperty<RealTensorValue>("secondPiolaKirchhofTensorBar")),
    _Siso(declareProperty<RealTensorValue>("secondPiolaKirchhofTensorIso")),
    _Pdev(declareProperty<RealTensorValue>("Pdev"))
   // _function2(getFunction("function"))
    
{
}

void
IncompressibleNeoHookean::computeQpPropertiesDerived()
{
  // here you have available _F, _U, _J because GuccioneCosta derives
  // from ElasticMaterial

  RealTensorValue const &F = _F[_qp];

  if (_mesh.dimension() == 3)
  {
    _C[_qp] = F.transpose() * F;
    _invC[_qp] = _C[_qp].inverse();
  }
  if (_mesh.dimension() == 2)
  {
    RealTensorValue Ftemp = F;
    Ftemp(2, 2) = 1.0;

    _C[_qp] = Ftemp.transpose() * Ftemp;
    _invC[_qp] = _C[_qp].inverse();

    for (int i = 0; i < 3; ++i)
    {

      _C[_qp](i, 2) = 0.0;
      _C[_qp](2, i) = 0.0;

      _invC[_qp](i, 2) = 0.0;
      _invC[_qp](2, i) = 0.0;
    }
  }
  if (_mesh.dimension() == 1)
  {
    _C[_qp] = F(0, 0) * F(0, 0);
    _invC[_qp](0, 0) = 1. / _C[_qp](0, 0);
  }

  _J23[_qp] = std::pow(_J[_qp], -2. / _mesh.dimension());

  _Sbar[_qp] = mu * _Id;

  _Siso[_qp] = _J23[_qp] * (_Sbar[_qp] -
                            (1.0 / _mesh.dimension()) *
                                _Sbar[_qp].contract(_C[_qp]) * _invC[_qp]);

  _P[_qp] += F * _Siso[_qp];
    
  _Pdev[_qp]=_P[_qp];

}

void
IncompressibleNeoHookean::evaluate_stress_lin(unsigned const &qp,
                                              RealTensorValue const &H,
                                              RealTensorValue &stressLin)
{

  RealTensorValue const &F = _F[qp];
  Clin = (F.transpose() * H + H.transpose() * F);

  stressLin =
      -1.0 / _mesh.dimension()  * _invC[qp].contract(Clin) * _Pdev[qp] 
      +  H * _Siso[qp] + _J23[qp] * _F[qp] * (-1.0 / _mesh.dimension() * _Sbar[qp].contract(Clin) * _invC[qp] +
                           1.0 / _mesh.dimension() * _Sbar[qp].contract(_C[qp]) * _invC[qp] *
                               Clin * _invC[qp]);
}
