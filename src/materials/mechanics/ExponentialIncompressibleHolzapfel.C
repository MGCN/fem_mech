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
#include "ExponentialIncompressibleHolzapfel.h"
#include "MooseMesh.h"

registerMooseObject("MECHApp", ExponentialIncompressibleHolzapfel);

template <>
InputParameters
validParams<ExponentialIncompressibleHolzapfel>()
{
  InputParameters params = validParams<ElasticMaterial>();
  params.addRequiredParam<Real>("c10", "c10");
  params.addRequiredParam<Real>("c01", "c01");
  params.addRequiredParam<Real>("k_1f", "k_1f");
  params.addRequiredParam<Real>("k_2f", "k_2f");
  params.addRequiredParam<Real>("kappa", "kappa");
  return params;
}

ExponentialIncompressibleHolzapfel::ExponentialIncompressibleHolzapfel(InputParameters const &parameters)
    : ElasticMaterial(parameters), 
      c01(getParam<Real>("c01")),
      c10(getParam<Real>("c10")),
      k_1f(getParam<Real>("k_1f")),
      k_2f(getParam<Real>("k_2f")),


      _I1(declareProperty<Real>("firstInvariant")),
      _I4(declareProperty<Real>("fourthInvariant")),
      _I1_bar(declareProperty<Real>("firstInavriantbar")),
      _I4_bar(declareProperty<Real>("fourthInvariantfbar")),
     
      _c1(declareProperty<Real>("coeffs1")),
      _c4_iso(declareProperty<Real>("coeffs2")),
      _c4_aniso(declareProperty<Real>("coeffs2")),
      kappa(declareProperty<Real>("kappa")),
  

      _fXf(getMaterialProperty<RealTensorValue>("fiberOuterFiber")),

      _C(declareProperty<RealTensorValue>("CauchyGreen")),
      _invC(declareProperty<RealTensorValue>("inverseCauchyGreen")),
      _J23(declareProperty<Real>("J23")),
      
      _S1bar(declareProperty<RealTensorValue>("_S1bar")),
      _S4bar(declareProperty<RealTensorValue>("_S4bar")),
      _Sbar(declareProperty<RealTensorValue>("secondPiolaKirchhofTensorBar")),
      _Siso(declareProperty<RealTensorValue>("secondPiolaKirchhofTensorIso")),
      _Pdev(declareProperty<RealTensorValue>("Pdev"))


{
}

void
ExponentialIncompressibleHolzapfel::computeQpPropertiesDerived()
{
  
  RealTensorValue const &F = _F[_qp];
   
  _C[_qp] = F.transpose() * F;
  
  RealTensorValue const &C = _C[_qp];
 
  _invC[_qp] = _C[_qp].inverse();

  _J23[_qp] = std::pow(_J[_qp], -2. / _mesh.dimension());

  _I1[_qp] = C.tr();

  _I4[_qp] = C.contract(_fXf[_qp]);
  
  _I1_bar[_qp] = _J23[_qp] * _I1[_qp];

  _I4_bar[_qp] = _J23[_qp] * _I4[_qp];
  
  Real const &I4 = _I4[_qp];
  
  Real const &I4_bar = _I4_bar[_qp];

  Real const &I1_bar = _I1_bar[_qp];


  //Isotropic Contribution

  _c1[_qp] = c10 * c01 * std::exp(c01 * (I1_bar - 3.0));

  Real const &c1 = _c1[_qp];
  
  //Anisotropic Contribution

  //N.B the fourth invariant acts just if you have a tensile state

  if (I4_bar > 1.0)
    _c4_aniso[_qp] = 2.0 * k_1f  * ( kappa[_qp] * I1_bar + (1.0 - 3.0 * kappa[_qp]) * I4_bar - 1.0 ) * ( std::exp(k_2f * std::pow(kappa[_qp] * I1_bar + 
                     ( 1.0 - 3.0 * kappa[_qp]) * I4_bar - 1.0, 2.0))) * (1.0 - 3.0 * kappa[_qp]);

  else
    _c4_aniso[_qp] = 0.0;



  if (I4_bar > 1.0)
    _c4_iso[_qp] = 2.0 * k_1f  * (kappa[_qp] * I1_bar *  (1.0 - 3.0 * kappa[_qp]) * I4_bar - 1.0) * (std::exp(k_2f * std::pow(kappa[_qp] * I1_bar + 
                      (1.0 - 3.0 * kappa[_qp]) * I4_bar - 1.0, 2.0))) * kappa[_qp];
  else
    _c4_iso[_qp] = 0.0;                  


  Real const &c4_aniso = _c4_aniso[_qp];


  Real const &c4_iso = _c4_iso[_qp];


    
  _S1bar[_qp]  = c1 * _Id;

  _S4bar[_qp]  = c4_aniso * _fXf[_qp] +  c4_iso * _Id;

    
//  std::cout<<" _fXf[_qp]"<< _fXf[_qp]<<std::endl;
   
  _Sbar[_qp] = _S1bar[_qp] +  _S4bar[_qp];
  

  _Siso[_qp] = _J23[_qp] * (_Sbar[_qp] - (1.0 / _mesh.dimension()) * _Sbar[_qp].contract(_C[_qp]) * _invC[_qp]);
  
  
  _P[_qp] += F * _Siso[_qp];
    

  _Pdev[_qp] = _P[_qp];
}

void
ExponentialIncompressibleHolzapfel::evaluate_stress_lin(unsigned const &qp, RealTensorValue const &H,
                                   RealTensorValue &stressLin)
{

    compute_C_lin(qp, H);



    RealTensorValue _S1bar_lin(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), 
                    _S41_bar_lin(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                    _S44_bar_lin(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                    _S41_bar_lin_temp(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                    _S44_bar_lin_temp(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  

  Real const &I4 = _I4[_qp];
  
  Real const &I4_bar = _I4_bar[_qp];

  Real const &I1_bar = _I1_bar[_qp];

  Real const &I1 = _I1[_qp];
    
 // RealTensorValue temp = c10 * c01 * c01 * std::exp(c01 * (I1_bar - 3.0)) * _I1_bar_lin;  

 // _S1bar_lin = temp.contract(_CLin);


  RealTensorValue I1_bar_lin = -1.0/_mesh.dimension() * _J23[qp] *  _invC[qp] * _I1[qp] +  _J23[qp] * _Id;

  RealTensorValue I4_bar_lin = -1.0/_mesh.dimension() * _J23[qp] *  _invC[qp] * _I4[qp] +  _J23[qp] * _fXf[qp];
    
  RealTensorValue temp = c10 * c01 * c01 * std::exp(c01 * (I1_bar - 3.0)) * I1_bar_lin;

  _S1bar_lin = temp.contract(_CLin);   


  Cbar_lin = _J23[qp] *  ( - 1.0 / _mesh.dimension() * _invC[qp].contract(_CLin) * _C[qp] + _CLin);

  double expo= std::exp(k_2f * std::pow(kappa[_qp] * I1_bar + (1.0 - 3.0 * kappa[_qp]) * I4_bar - 1.0, 2.0));
          

  if (_I4_bar[qp] > 1.0)
  {
   /*  _S44_bar_lin_temp = ( 2.0 * k_1f * expo * ( 1 - 3 * kappa[_qp]) + 2 * k_2f * _c4_aniso[_qp] *  ( kappa[_qp] * I1_bar +
                        ( 1.0 - 3.0 * kappa[_qp]) * I4_bar - 1.0) ) * ( kappa[_qp] * I1_bar_lin + (1.0 - 3.0 * kappa[_qp]) * I4_bar_lin ); 


    _S41_bar_lin_temp = ( 2.0 * k_1f * expo *  kappa[_qp] + 2 * k_2f * _c4_iso[_qp] *  ( kappa[_qp] * I1_bar +
                        ( 1.0 - 3.0 * kappa[_qp]) * I4_bar - 1.0) ) * ( kappa[_qp] * I1_bar_lin + (1.0 - 3.0 * kappa[_qp]) * I4_bar_lin );*/
    _S44_bar_lin_temp = (2.0 * k_1f  * ( kappa[_qp] * I1_bar_lin + (1.0 - 3.0 * kappa[_qp]) * I4_bar_lin ) * ( std::exp(k_2f * std::pow(kappa[_qp] * I1_bar +
                        (1.0 - 3.0 * kappa[_qp]) * I4_bar - 1.0, 2.0))) * (1.0 - 3.0 * kappa[_qp]) + _c4_iso[qp] * k_2f * ( kappa[_qp] * I1_bar +
                        (1.0 - 3.0 * kappa[_qp]) * I4_bar - 1.0) * (kappa[_qp] * I1_bar_lin + (1-3.0 * kappa[_qp]) * I4_bar_lin) ) * kappa[_qp] *_Id +
                        (2.0 * k_1f  * ( kappa[_qp] * I1_bar_lin + (1.0 - 3.0 * kappa[_qp]) * I4_bar_lin ) * ( std::exp(k_2f * std::pow(kappa[_qp] * I1_bar +
                        (1.0 - 3.0 * kappa[_qp]) * I4_bar - 1.0, 2.0))) * (1.0 - 3.0 * kappa[_qp]) + _c4_aniso[qp] * 2 * k_2f * ( kappa[_qp] * I1_bar +
                        (1.0 - 3.0 * kappa[_qp]) * I4_bar - 1.0) * (kappa[_qp] * I1_bar_lin + (1- 3.0 * kappa[_qp]) * I4_bar_lin) * ( 1 - 3.0 * kappa[_qp])) * _fXf[qp];
 } 
 


  _S44_bar_lin = _S44_bar_lin_temp.contract(_CLin); 

  //_S41_bar_lin = _S41_bar_lin_temp.contract(_CLin); 
   
  _Sbar_lin =  _S1bar_lin +  _S44_bar_lin;

  stressLin =  -  1.0 / _mesh.dimension() *  _invC[qp].contract(_CLin) * _Pdev[qp] +  H * _Siso[qp] +
  
                +  _J23[qp] * 1.0 / _mesh.dimension() * _F[qp] *  (_mesh.dimension() * _Sbar_lin - _Sbar[qp].contract(_CLin) * _invC[qp] - _Sbar_lin.contract(_C[qp]) * _invC[qp] +
                                            
                                                                   + 1.0 * _Sbar[qp].contract(_C[qp]) * (_invC[qp] * _CLin * _invC[qp])); //(_invC[qp] * _CLin * _invC[qp]));
    
}

void
ExponentialIncompressibleHolzapfel::compute_C_lin(unsigned const &qp, RealTensorValue const &H)
{
   _CLin = (_F[qp].transpose() * H + H.transpose() * _F[qp]);
}

 
