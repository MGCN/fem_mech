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
#include "IncompressibleHolzapfel.h"
#include "MooseMesh.h"

registerMooseObject("MECHApp", IncompressibleHolzapfel);

template <>
InputParameters
validParams<IncompressibleHolzapfel>()
{
  InputParameters params = validParams<ElasticMaterial>();
  params.addRequiredParam<Real>("mu", "mu");
  params.addRequiredParam<Real>("k_1f", "k_1f");
  params.addRequiredParam<Real>("k_2f", "k_2f");
  params.addRequiredParam<Real>("k_1s", "k_1s");
  params.addRequiredParam<Real>("k_2s", "k_2s");
  params.addRequiredParam<Real>("k_1fs", "k_1fs");
  params.addRequiredParam<Real>("k_2fs", "k_2fs");
  return params;
}

IncompressibleHolzapfel::IncompressibleHolzapfel(InputParameters const &parameters)
    : ElasticMaterial(parameters), 
      mu(getParam<Real>("mu")),
      k_1f(getParam<Real>("k_1f")),
      k_2f(getParam<Real>("k_2f")),
      k_1s(getParam<Real>("k_1s")),
      k_2s(getParam<Real>("k_2s")),
      k_1fs(getParam<Real>("k_1fs")),
      k_2fs(getParam<Real>("k_2fs")),

      _I4(declareProperty<Real>("fourthInvariantf")),
      _I6(declareProperty<Real>("sixthInvariants")),
      _I8(declareProperty<Real>("eightInvariants")),
      _I4_bar(declareProperty<Real>("fourthInvariantfbar")),
      _I6_bar(declareProperty<Real>("sixthInvariantsbar")),
      _I8_bar(declareProperty<Real>("eightInvariantsbar")),
     
      _c1(declareProperty<Real>("coeffs1")),
      _c4(declareProperty<Real>("coeffs2")),
      _c6(declareProperty<Real>("coeffs3")),
      _c8(declareProperty<Real>("coeffs4")),
  

      _fXf(getMaterialProperty<RealTensorValue>("fiberOuterFiber")),
      _gXg(getMaterialProperty<RealTensorValue>("gfiberOutergFiber")), //gfiberOutergFiber
      _sXs(getMaterialProperty<RealTensorValue>("sheetOuterSheet")),
      _fXs(getMaterialProperty<RealTensorValue>("fiberOuterSheet")),

      _C(declareProperty<RealTensorValue>("CauchyGreen")),
      _invC(declareProperty<RealTensorValue>("inverseCauchyGreen")),
      _J23(declareProperty<Real>("J23")),
      
      _S1bar(declareProperty<RealTensorValue>("_S1bar")),
      _S4bar(declareProperty<RealTensorValue>("_S2bar")),
      _S6bar(declareProperty<RealTensorValue>("_S3bar")),
      _S8bar(declareProperty<RealTensorValue>("_S3bar")),
      _Sbar(declareProperty<RealTensorValue>("secondPiolaKirchhofTensorBar")),
      _Siso(declareProperty<RealTensorValue>("secondPiolaKirchhofTensorIso")),
      _Pdev(declareProperty<RealTensorValue>("Pdev"))


{
}

void
IncompressibleHolzapfel::computeQpPropertiesDerived()
{
  
  RealTensorValue const &F = _F[_qp];
   
  _C[_qp] = F.transpose() * F;
  
  RealTensorValue const &C = _C[_qp];
 
  _invC[_qp] = _C[_qp].inverse();

  _J23[_qp] = std::pow(_J[_qp], -2. / _mesh.dimension());

  _I4[_qp] = C.contract(_fXf[_qp]);
  
  _I4_bar[_qp] = _J23[_qp] * _I4[_qp];
  
  Real const &I4 = _I4[_qp];
  
  Real const &I4_bar = _I4_bar[_qp];

  _I6[_qp] = C.contract(_gXg[_qp]);
  
  _I6_bar[_qp] = _J23[_qp] * _I6[_qp];
  
  Real const &I6 = _I6[_qp];
  
  Real const &I6_bar = _I6_bar[_qp];
    
  _I8[_qp] = C.contract(_fXs[_qp]);
    
  _I8_bar[_qp] = _J23[_qp] * _I8[_qp];
    
  Real const &I8 = _I8[_qp];
    
  Real const &I8_bar = _I8_bar[_qp];

  //Isotropic Contribution

  _c1[_qp] = mu ;

  Real const &c1 = _c1[_qp];
  
  //Anisotropic Contribution

  //N.B the fourth invariant acts just if you have a tensile state

  if (I4_bar > 1.0)
    _c4[_qp] =
        2.0 * k_1f * (std::exp(k_2f * std::pow(_I4_bar[_qp] - 1.0, 2.0)) * (I4_bar - 1.0));

  else
    _c4[_qp] = 0.0;

  Real const &c4 = _c4[_qp];


  if (I6_bar > 1.0)

    _c6[_qp] =
        2.0 * k_1s * (std::exp(k_2s * std::pow(_I6_bar[_qp] - 1.0, 2.0)) * (I6_bar - 1.0));
  else

    _c6[_qp] = 0.0;

  Real const &c6 = _c6[_qp];
    
  _c8[_qp] =
      2.0 * k_1fs * std::exp(k_2fs * I8_bar * I8_bar) * (I8_bar);

  Real const &c8 = _c8[_qp];
    
  _S1bar[_qp] = c1 * _Id;
  _S4bar[_qp] = c4 * _fXf[_qp];
  _S6bar[_qp] = c6 * _gXg[_qp];
  _S8bar[_qp] = c8 * _fXs[_qp];
    
//  std::cout<<" _fXf[_qp]"<< _fXf[_qp]<<std::endl;
   
  _Sbar[_qp] = _S1bar[_qp] +  _S4bar[_qp] +  _S6bar[_qp] +  _S8bar[_qp];
  

  _Siso[_qp] = _J23[_qp] * (_Sbar[_qp] - (1.0 / _mesh.dimension()) * _Sbar[_qp].contract(_C[_qp]) * _invC[_qp]);
  
  
  _P[_qp] += F * _Siso[_qp];
    

  _Pdev[_qp] = _P[_qp];
}

void
IncompressibleHolzapfel::evaluate_stress_lin(unsigned const &qp, RealTensorValue const &H,
                                   RealTensorValue &stressLin)
{

    compute_C_lin(qp, H);



    RealTensorValue _S1bar_lin(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), 
                    _S4bar_lin(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                    _S6bar_lin(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                    _S8bar_lin(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  
    
    _S1bar_lin = 0;  
    
   
    Cbar_lin = _J23[qp] *  ( - 1.0 / _mesh.dimension() * _invC[qp].contract(_CLin) * _C[qp] + _CLin);
    
         

  if (_I4_bar[qp] > 1.0)
  {
    _S4bar_lin = (2.0 * k_2f * _c4[qp] * (_I4_bar[qp] - 1.0) +
                  2.0 * k_1f * std::exp(k_2f * std::pow(_I4_bar[qp] - 1.0, 2.0))) *
                 _J23[qp] * (_fXf[qp].contract(_CLin) +  _I4[qp] * ( - 1.0 / _mesh.dimension() * _invC[qp].contract(_CLin))) * _fXf[qp];
  } 
  
  

  if (_I6_bar[qp] > 1.0)
  {
    _S6bar_lin = (2.0 * k_2s * _c6[qp] * (_I6_bar[qp] - 1.0) +
                  2.0 * k_1s * std::exp(k_2s * std::pow(_I6_bar[qp] - 1.0, 2.0)))*
                 _J23[qp] * (_gXg[qp].contract(_CLin) + _I6[qp] * ( - 1.0 / _mesh.dimension() * _invC[qp].contract(_CLin))) * _gXg[qp];
      
  }

    

    _S8bar_lin = (2.0 * k_2fs * _c8[qp] * (_I8_bar[qp]) +
                  2.0 * k_1fs * std::exp(k_2fs * std::pow(_I8_bar[qp],2.0)))*
                  _J23[qp] * (_fXs[qp].contract(_CLin) + _I8[qp] * ( - 1.0 / _mesh.dimension() * _invC[qp].contract(_CLin)))* _fXs[qp];
    



     
   
    _Sbar_lin =  _S1bar_lin + _S4bar_lin + _S6bar_lin + _S8bar_lin;

    stressLin =  -  1.0 / _mesh.dimension() *  _invC[qp].contract(_CLin) * _Pdev[qp] +  H * _Siso[qp] +
      
                +  _J23[qp] * 1.0 / _mesh.dimension() * _F[qp] *  (_mesh.dimension() * _Sbar_lin - _Sbar[qp].contract(_CLin) * _invC[qp] - _Sbar_lin.contract(_C[qp]) * _invC[qp] +
                                            
                                                                   + 1.0 * _Sbar[qp].contract(_C[qp]) * (_invC[qp] * _CLin * _invC[qp])); //(_invC[qp] * _CLin * _invC[qp]));
    
}

void
IncompressibleHolzapfel::compute_C_lin(unsigned const &qp, RealTensorValue const &H)
{
   _CLin = (_F[qp].transpose() * H + H.transpose() * _F[qp]);
}

 
