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
#include "IncompressibleMooneyRivlin.h"
#include "MooseMesh.h"

registerMooseObject("MECHApp", IncompressibleMooneyRivlin);

template <>
InputParameters
validParams<IncompressibleMooneyRivlin>()
{
  InputParameters params = validParams<ElasticMaterial>();
  params.addRequiredParam<Real>("c10", "c10");
  params.addRequiredParam<Real>("c20", "c20");
  return params;
}

IncompressibleMooneyRivlin::IncompressibleMooneyRivlin(
    const InputParameters &parameters)
    : ElasticMaterial(parameters), 
      _c10(getParam<Real>("c10")),
      _c20(getParam<Real>("c20")),
      _C(declareProperty<RealTensorValue>("CauchyGreen")),
      _invC(declareProperty<RealTensorValue>("inverseCauchyGreen")),
      _J23(declareProperty<Real>("J23")),
      _Sbar(declareProperty<RealTensorValue>("secondPiolaKirchhofTensorBar")),
      _Siso(declareProperty<RealTensorValue>("secondPiolaKirchhofTensorIso"))
{
}

void
IncompressibleMooneyRivlin::computeQpPropertiesDerived()
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

    _Sbar[_qp] = _c10 * _Id + _c20 * ((_J23[_qp] * _C[_qp]).tr() * _Id - _J23[_qp]*  _C[_qp] );

    _Siso[_qp] = _J23[_qp] * (_Sbar[_qp] - (1.0 / _mesh.dimension()) * _Sbar[_qp].contract(_C[_qp]) * _invC[_qp]);

    _P[_qp] += F * _Siso[_qp];
}

void
IncompressibleMooneyRivlin::evaluate_stress_lin(unsigned const &qp,
                                              RealTensorValue const &H,
                                              RealTensorValue &stressLin)
{

   RealTensorValue const &F = _F[qp];
  
   Clin = (F.transpose() * H + H.transpose() * F);
   
   elin = ( H + H.transpose());
  
   RealTensorValue  C_bar_lin = _J23[qp] *  ( - 1.0 / _mesh.dimension() * _invC[qp].contract(Clin) * _C[qp] + Clin);
     
   _Sbar_lin = _c20 * _Id.contract(C_bar_lin) * _Id - _c20 * C_bar_lin ;

   //_Sbar_lin = - 1.0/3.0 * _J23[qp] * _c20 * _invC[qp].contract(Clin) * _Id  + _c20 * _J23 [_qp] * _Id.contract(Clin) * _Id - _c20 * _J23[qp] * (Clin - 1.0 / 3.0 * _invC[qp].contract(Clin) * _C[qp]);
   

   _Siso_lin =  - _J23[qp] * 1.0 / _mesh.dimension() *  _invC[qp].contract(Clin) * (_Sbar[qp] - (1.0 / _mesh.dimension()) * _Sbar[qp].contract(_C[qp]) * _invC[qp]) +
  
                 _J23[qp] * 1.0 / _mesh.dimension() *  ( _mesh.dimension() * _Sbar_lin - _Sbar[qp].contract(Clin) * _invC[qp] - _Sbar_lin.contract(_C[qp]) * _invC[qp] +
                                            
                                           _Sbar[qp].contract(_C[qp]) * _invC[qp] * Clin * _invC[qp]);                                     
                                        
  stressLin = F * _Siso_lin + H * _Siso[qp];  
                                           
}

