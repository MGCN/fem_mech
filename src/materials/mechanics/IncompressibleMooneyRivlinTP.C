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
#include "IncompressibleMooneyRivlinTP.h"
#include "MooseMesh.h"

registerMooseObject("MECHApp", IncompressibleMooneyRivlinTP);

template <>
InputParameters
validParams<IncompressibleMooneyRivlinTP>()
{
    InputParameters params = validParams<ElasticMaterial>();
    params.addRequiredParam<Real>("c10", "c10");
    params.addRequiredParam<Real>("c20", "c20");
    params.addRequiredParam<Real>("c11", "c11");
    
    return params;
}

IncompressibleMooneyRivlinTP::IncompressibleMooneyRivlinTP(
                           const InputParameters &parameters)
: ElasticMaterial(parameters),
  _c10(getParam<Real>("c10")),
  _c20(getParam<Real>("c20")),
  _c11(getParam<Real>("c11")),
  _I1_bar(declareProperty<Real>("firstInvariant")),
  _I2_bar(declareProperty<Real>("secondInvariant")),
  _C(declareProperty<RealTensorValue>("CauchyGreen")),
  _C_bar(declareProperty<RealTensorValue>("CauchyGreenBar")),
  _invC(declareProperty<RealTensorValue>("inverseCauchyGreen")),
  _J23(declareProperty<Real>("J23")),
  _Sbar(declareProperty<RealTensorValue>("secondPiolaKirchhofTensorBar")),
  _Siso(declareProperty<RealTensorValue>("secondPiolaKirchhofTensorIso"))
{

}



void
IncompressibleMooneyRivlinTP::computeQpPropertiesDerived()
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


    _C_bar[_qp] = _J23[_qp] * _C[_qp];

    RealTensorValue const &C_bar_sq = _C_bar[_qp] * _C_bar[_qp];

    _I1_bar[_qp] = _C_bar[_qp].tr();

    _I2_bar[_qp]= 1.0/2.0 * (_C[_qp].tr() * _C[_qp].tr() - C_bar_sq.tr()) ;

    _Sbar[_qp] = _c10 * _Id + _c20 * (_I1_bar[_qp] * _Id -  _C_bar[_qp]) + _c11 * (_I2_bar[_qp] - 1.0) * _Id + _c11 * (_I1_bar[_qp] - 1.0) * ( _I1_bar[_qp] * _Id  -  _C_bar[_qp]);

    _Siso[_qp] = 2 * _J23[_qp] * (_Sbar[_qp] - (1.0 / _mesh.dimension()) * _Sbar[_qp].contract(_C[_qp]) * _invC[_qp]);

    _P[_qp] += F * _Siso[_qp];
}

void
IncompressibleMooneyRivlinTP::evaluate_stress_lin(unsigned const &qp,
                                              RealTensorValue const &H,
                                              RealTensorValue &stressLin)
{

  //  RealTensorValue const &F = _F[qp];
  
   
  //  Clin = (F.transpose() * H + H.transpose() * F);

  
  //  RealTensorValue  C_bar_lin =  _J23[qp] *  ( - 1.0 / _mesh.dimension() * _invC[qp].contract(Clin) * _C[qp] + Clin);


  //  RealTensorValue const &C_bar = _J23[qp] * _C[qp];


  //  Real const &I_1_bar = C_bar.tr();
     
   
  //  RealTensorValue  const &_term_1 = _Id.contract(C_bar_lin) * _Id - C_bar_lin ;

  
  //  RealTensorValue  const &_term_2 =  (_I1_bar[qp] * _Id -  _C_bar[qp]).contract(C_bar_lin) * _Id;
     
  
  //  _Sbar_lin = _c20 * _term_1 + _c11 * _term_2 + _c11 * _term_1 * (_I1_bar[qp] - 1.0) + _c11 * (_I1_bar[qp] * _Id  - _C_bar[qp]) * _Id.contract(C_bar_lin);

  //  // + _c11 * _term_3;

  //  //_Id.contract(C_bar_lin) * _Id - _c20 * C_bar_lin ;
  

  //  _Siso_lin =  - _J23[qp] * 1.0 / _mesh.dimension() *  _invC[qp].contract(Clin) * (_Sbar[qp] - (1.0 / _mesh.dimension()) * _Sbar[qp].contract(_C[qp]) * _invC[qp]) +
  
  //                _J23[qp] * 1.0 / _mesh.dimension() *  ( _mesh.dimension() * _Sbar_lin - _Sbar[qp].contract(Clin) * _invC[qp] - _Sbar_lin.contract(_C[qp]) * _invC[qp] +
                                            
  //                                          _Sbar[qp].contract(_C[qp]) * _invC[qp] * Clin * _invC[qp]);                                     
                                                                          
                                        
  // stressLin = 2 * F * _Siso_lin + 2 * H * _Siso[qp];  

  RealTensorValue const &F = _F[qp];
  
  Clin = (F.transpose() * H + H.transpose() * F);
   
  elin = ( H + H.transpose());
  
  RealTensorValue  C_bar_lin = _J23[qp] *  ( - 1.0 / _mesh.dimension() * _invC[qp].contract(Clin) * _C[qp] + Clin);
     
  RealTensorValue  const & _Sbar_lin_c_20 = _c20 * _Id.contract(C_bar_lin) * _Id - _c20 * C_bar_lin ;


  RealTensorValue  const &  _Sbar_lin_c_11 = _c11 * (_Id.contract(C_bar_lin) * _Id - C_bar_lin) * (_I1_bar[qp] - 1.0) + _c11 *( _I1_bar[qp] * _Id -  _C_bar[qp]) * _Id.contract(C_bar_lin) + 

  _c11 * (_I1_bar[qp] * _Id - _C_bar[qp]).contract(C_bar_lin) * _Id;


  _Sbar_lin = _Sbar_lin_c_20 + _Sbar_lin_c_11;

   //_Sbar_lin = - 1.0/3.0 * _J23[qp] * _c20 * _invC[qp].contract(Clin) * _Id  + _c20 * _J23 [_qp] * _Id.contract(Clin) * _Id - _c20 * _J23[qp] * (Clin - 1.0 / 3.0 * _invC[qp].contract(Clin) * _C[qp]);
   

   _Siso_lin =  - _J23[qp] * 1.0 / _mesh.dimension() *  _invC[qp].contract(Clin) * (_Sbar[qp] - (1.0 / _mesh.dimension()) * _Sbar[qp].contract(_C[qp]) * _invC[qp]) +
  
                 _J23[qp] * 1.0 / _mesh.dimension() *  ( _mesh.dimension() * _Sbar_lin - _Sbar[qp].contract(Clin) * _invC[qp] - _Sbar_lin.contract(_C[qp]) * _invC[qp] +
                                            
                                           _Sbar[qp].contract(_C[qp]) * _invC[qp] * Clin * _invC[qp]);                                     
                                        
  stressLin = 2 * F * _Siso_lin + 2 * H * _Siso[qp]; 
                                           
}




