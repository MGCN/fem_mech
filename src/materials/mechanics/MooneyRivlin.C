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
/*                Prepared by Maria Nestola,                    */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/****************************************************************/
#include "MooneyRivlin.h"
#include "MooseMesh.h"

registerMooseObject("MECHApp", MooneyRivlin);

template <>
InputParameters
validParams<MooneyRivlin>()
{
    InputParameters params = validParams<ElasticMaterial>();
    params.addRequiredParam<Real>("c10", "c10");
    params.addRequiredParam<Real>("c20", "c20");
    params.addRequiredParam<Real>("lambda", "lambda");
    
    return params;
}

MooneyRivlin::MooneyRivlin(
                           const InputParameters &parameters)
: ElasticMaterial(parameters),
_c10(getParam<Real>("c10")),
_c20(getParam<Real>("c20")),
_lambda(getParam<Real>("lambda")),
_C(declareProperty<RealTensorValue>("CauchyGreen")),
_invC(declareProperty<RealTensorValue>("inverseCauchyGreen"))
{
}

void
MooneyRivlin::computeQpPropertiesDerived()
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
    
    
    //_P[_qp] +=  1.0 * _c10 * (F - _invFtr[_qp]) + _lambda * std::log(_J[_qp]) * _invFtr[_qp];
    
    _P[_qp] +=  2.0 * _c10 * F + 2.0 * _c20 * F * ( _C[_qp] * _Id - _C[_qp] )  - (2.0 * _c10 + 4.0 * _c20) *  _invFtr[_qp] + _lambda * std::log(_J[_qp]) * _invFtr[_qp];
}

void
MooneyRivlin::evaluate_stress_lin(unsigned const &qp,
                                  RealTensorValue const &H,
                                  RealTensorValue &stressLin)
{
    
    Real const &J = _J[qp];
    
    RealTensorValue const &invFtr = _invFtr[qp];
    
    RealTensorValue const &F = _F[_qp];
    
    //RealTensorValue Plin = _c10 * H - (1.0 * _lambda * std::log(J) - 1.0 * _c10) * invFtr * H.transpose() * invFtr + _lambda * invFtr.contract(H) * invFtr;
    
    RealTensorValue Plin = 2.0 * _c10 * H + 2.0 * _c20 * H * (  _C[_qp] * _Id - _C[_qp] ) - (1.0 * _lambda * std::log(J) - 1.0 * (2.0 * _c10 + 4.0 * _c20)) * invFtr * H.transpose() * invFtr
                           + _lambda * invFtr.contract(H) * invFtr;
    
    stressLin = Plin;
    
}

