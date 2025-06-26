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
/*               Mech  - ICS Mechanical simulation framework    */
/*                Prepared by Maria Nestola                     */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/*.                                                             */
/****************************************************************/
#include "GrowthMaterial.h"
registerMooseObject("MECHApp", GrowthMaterial);

template <>

InputParameters
validParams<GrowthMaterial>()
{
  InputParameters params = validParams<ElasticMaterial>();
  params.addRequiredParam<Real>("mu", "mu");
  params.addRequiredParam<Real>("lambda", "lambda");
  params.addRequiredCoupledVar("theta","theta");
  return params;
}

GrowthMaterial::GrowthMaterial(const InputParameters &parameters)
    : ElasticMaterial(parameters), 
      _mu(getParam<Real>("mu")),
      _lambda(getParam<Real>("lambda")),
      _theta(coupledValue("theta")),
      _C(declareProperty<RealTensorValue>("CauchyGreen")),
      _invC(declareProperty<RealTensorValue>("inverseCauchyGreen")),
      _S(declareProperty<RealTensorValue>("SecondPiola")),
      // _Ce(declareProperty<RealTensorValue>("CauchyGreenElastic")),
      _invCe(declareProperty<RealTensorValue>("InverseCauchyGreenElastic")),
      // _Se(declareProperty<RealTensorValue>("SecondPiolaElastic")),
      _Fg(declareProperty<RealTensorValue>("GrowthTensor")),
      _Fgt(declareProperty<RealTensorValue>("TransposeGrowthTensor")),
      _invFg(declareProperty<RealTensorValue>("inverseGrowthTensor")),
      _invFgt(declareProperty<RealTensorValue>("InverseTransposeGrowthTensor"))
{
}

void
GrowthMaterial::computeQpPropertiesDerived()
{
  // here you have available _F, _U, _stress, J because NeoHookean derives
  // from ElasticMaterial

    RealTensorValue const &F = _F[_qp];



    RealTensorValue const &Fg = _theta[_qp] * _Id;
    
    if (_mesh.dimension() == 3)
    {
      _C[_qp]      = F.transpose() * F;

      _invC[_qp]   = _C[_qp].inverse();

      _Fgt[_qp]    = Fg.transpose();

      _invFg[_qp]  = Fg.inverse();

      _invFgt[_qp] = _Fgt[_qp].inverse();

      _Ce[_qp]     = _invFgt[_qp] * _C[_qp] * _invFg[_qp];

      _invCe[_qp]  = _Ce[_qp].inverse();
    }
    
    if (_mesh.dimension() == 2)
    {
        RealTensorValue Ftemp = F;

        RealTensorValue Fgtemp = Fg;
        
        Ftemp(2, 2)  = 1.0;

        Fgtemp(2, 2) = 1.0;
        
        _C[_qp]      = Ftemp.transpose() * Ftemp;

        _invC[_qp]   = _C[_qp].inverse();

        _Fgt[_qp]    = Fgtemp.transpose();

        _invFg[_qp]  = Fgtemp.inverse();

        _invFgt[_qp] = _Fgt[_qp].inverse();

        _Ce[_qp]        = _invFgt[_qp] * _C[_qp] * _invFg[_qp];

        _invCe[_qp]     = _Ce[_qp].inverse();
        
        for (int i = 0; i < 3; ++i)
        {
            
            _C[_qp](i, 2) = 0.0;
            _C[_qp](2, i) = 0.0;
            
            _invC[_qp](i, 2) = 0.0;
            _invC[_qp](2, i) = 0.0;

            _Fg[_qp](i, 2) = 0.0;
            _Fg[_qp](2, i) = 0.0;
            
            _Ce[_qp](i, 2) = 0.0;
            _Ce[_qp](2, i) = 0.0;

            _invCe[_qp](i, 2) = 0.0;
            _invCe[_qp](2, i) = 0.0;
            
            _invFg[_qp](i, 2) = 0.0;
            _invFg[_qp](2, i) = 0.0;

            _Fgt[_qp](i, 2) = 0.0;
            _Fgt[_qp](2, i) = 0.0;
            
            _invFgt[_qp](i, 2) = 0.0;
            _invFgt[_qp](2, i) = 0.0;
        }
    }
    
    if (_mesh.dimension() == 1)
    {
      _C[_qp] = F(0, 0) * F(0, 0);

      _invC[_qp](0, 0) = 1. / _C[_qp](0, 0);

      _Fg[_qp] = _theta[_qp];

      _Fgt[_qp] = _theta[_qp];

      _invFg[_qp] = 1.0/_theta[_qp];

      _invFgt[_qp] = 1.0/_theta[_qp];

      _Ce[_qp] = _invFgt[_qp](0,0) * _C[_qp](0,0) * _invFg[_qp](0,0);

      _invCe[_qp] = 1.0 /_Ce[_qp](0,0);


    }
    

    Real const &J = _J[_qp];

    
    Real const &Je = 1.0/_theta[_qp] * _J[_qp];


    Real const &I3e = Je * Je;


    _Se[_qp] = _mu * _Id + (1.0/2 * _lambda * std::log(_Ce[_qp].det()) - _mu) * _invCe[_qp];

    _S[_qp] =  _invFg[_qp] * _Se[_qp] * _invFgt[_qp];

    _P[_qp] = 1.0/_theta[_qp] * 1.0/_theta[_qp] * _mu * F + (1.0/2.0 * _lambda * std::log(J * J) - _mu) * _invFtr[_qp];

}

void
GrowthMaterial::evaluate_stress_lin(unsigned const &qp, RealTensorValue const &H,
                                RealTensorValue &stressLin)

{

  Real const &J = _J[qp];

  RealTensorValue const &invFtr = _invFtr[qp];

  Real const &Je = 1.0/_theta[qp] * _J[qp];

  Real const &I3e = Je * Je;

  // RealTensorValue const &F = _F[_qp];

  // RealTensorValue  _CeLin = /*1.0/_theta[qp] * 1.0/_theta[qp]*/ 1.0 * (_F[qp].transpose() * H + H.transpose() * _F[qp]);

  //RealTensorValue Plin = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  
  //RealTensorValue  temp =  1.0/2.0 * _lambda * _invCe[qp] * (_invCe[qp].contract(_CeLin)) + (1.0/2 * _lambda * std::log( _Ce[qp].det()) - _mu) * (-1.0 * _invCe[qp] * _CeLin * _invCe[qp]);

  //stressLin = F * temp * /*1.0/_theta[qp] * 1.0/_theta[qp]*/ +  _mu * H + (1.0/2 * _lambda * std::log( _Ce[qp].det()) - _mu) * _invCe[qp] * H;

  RealTensorValue Plin = 1.0/_theta[qp] * 1.0/_theta[qp] * _mu * H - (1.0/2.0 * _lambda * std::log(J * J) - 1.0 * _mu) * invFtr * H.transpose() * invFtr + _theta[qp] * _theta[qp] * _lambda * invFtr.contract(H) * invFtr;
 //Plin;

  stressLin = Plin;
}



