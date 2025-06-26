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
/*   Immersed_Boundary - ICS Mechanical simulation framework    */
/*                Prepared by Maria Nestola                     */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/* Material to compute fibers direction.                        */
/****************************************************************/

#include "CompressibleHolzapfel.h"
#include "MooseMesh.h"

registerMooseObject("MECHApp",CompressibleHolzapfel);

template <>
InputParameters
validParams<CompressibleHolzapfel>()
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

CompressibleHolzapfel::CompressibleHolzapfel(InputParameters const &parameters)
    : ElasticMaterial(parameters), 
      mu(getParam<Real>("mu")),
      k_1f(getParam<Real>("k_1f")),
      k_2f(getParam<Real>("k_2f")),
      k_1s(getParam<Real>("k_1s")),
      k_2s(getParam<Real>("k_2s")),
      k_1fs(getParam<Real>("k_1s")),
      k_2fs(getParam<Real>("k_2s")),

      _I4(declareProperty<Real>("fourthInvariantf")),
      _I6(declareProperty<Real>("sixthInvariants")),
      _I8(declareProperty<Real>("sixthInvariants")),

     
      _c1(declareProperty<Real>("coeffs1")),
      _c4(declareProperty<Real>("coeffs2")),
      _c6(declareProperty<Real>("coeffs3")),
      _c8(declareProperty<Real>("coeffs4")),
  

      _fXf(getMaterialProperty<RealTensorValue>("fiberOuterFiber")),
      _sXs(getMaterialProperty<RealTensorValue>("sheetOuterSheet")),
      _fXs(getMaterialProperty<RealTensorValue>("fiberOuterSheet")),

      _C(declareProperty<RealTensorValue>("CauchyGreen")),
      _invC(declareProperty<RealTensorValue>("inverseCauchyGreen")),

      
      _S1(declareProperty<RealTensorValue>("_S1bar")),
      _S4(declareProperty<RealTensorValue>("_S2bar")),
      _S6(declareProperty<RealTensorValue>("_S3bar")),
      _S8(declareProperty<RealTensorValue>("_S3bar")),
      _S(declareProperty<RealTensorValue>("secondPiolaKirchhofTensorBar")),
      _Pdev(declareProperty<RealTensorValue>("Pdev"))


{
}

void
CompressibleHolzapfel::computeQpPropertiesDerived()
{
  
  RealTensorValue const &F = _F[_qp];
   
  _C[_qp] = F.transpose() * F;
  
  RealTensorValue const &C = _C[_qp];
 
  _invC[_qp] = _C[_qp].inverse();

  _I4[_qp] = C.contract(_fXf[_qp]);
 
  Real const &I4 = _I4[_qp];

  _I6[_qp] = C.contract(_sXs[_qp]);
  
  Real const &I6 = _I6[_qp];
    
  _I8[_qp] = C.contract(_fXs[_qp]);

  Real const &I8 = _I8[_qp];


  //Isotropic Contribution


  _c1[_qp] = mu ;

  Real const &c1 = _c1[_qp];
  
  //Anisotropic Contribution

  //N.B the fourth invariant acts just if you have a tensile state

  if (I4 > 1.0)
    _c4[_qp] =
        4.0 * k_1f/(2 * k_2f) * (std::exp(k_2f * (I4 - 1.0) * (I4 - 1.0)) * (I4 - 1.0));

  else
    _c4[_qp] = 0.0;

  Real const &c4 = _c4[_qp];


  if (I6 > 1.0)

    _c6[_qp] =
        4.0 * k_1s/(2 * k_2s) * std::exp(k_2s * (I6 - 1.0) * (I6 - 1.0)) * (I6 - 1.0);
  else
    _c6[_qp] = 0.0;

  Real const &c6 = _c6[_qp];
    
    _c8[_qp] =
      4.0 * k_1fs/(2 * k_2fs) * std::exp(k_2fs * I8 * I8) * (I8);

  Real const &c8 = _c8[_qp];
    
  _S1[_qp] = c1 * _Id;
  _S4[_qp] = c4 * _fXf[_qp];
  _S6[_qp] = c6 * _sXs[_qp];
  _S8[_qp] = c8 * _fXs[_qp];
    
    

   
  _S[_qp] = _S1[_qp] +  _S4[_qp] +  _S6[_qp] +  _S8[_qp];
  
  
  _P[_qp] += F * _S[_qp];

  _Pdev[_qp]=_P[_qp];
}

void
CompressibleHolzapfel::evaluate_stress_lin(unsigned const &qp, RealTensorValue const &H,
                                   RealTensorValue &stressLin)
{

    compute_C_lin(qp, H);



    RealTensorValue _S1_lin(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                    _S4_lin(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                    _S6_lin(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                    _S8_lin(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  
    
    _S1_lin = 0;
    
    
         

  if (_I4[qp] > 1.0)
  {
    _S4_lin = (2.0 * k_2f * _c4[qp] * (_I4[qp] - 1.0) +
                  4.0 * k_1f/(2 * k_2f) * std::exp(k_2f * std::pow(_I4[qp] - 1.0, 2.0))) * _fXf[qp].contract(_CLin) * _fXf[qp];
  } 
  
  

  if (_I6[qp] > 1.0)
  {
    _S6_lin = (2.0 * k_2s * _c6[qp] * (_I6[qp] - 1.0) +
                  4.0 * k_1s/(2 * k_2s) * std::exp(k_2s * std::pow(_I6[qp] - 1.0, 2.0))) * _sXs[qp].contract(_CLin) * _sXs[qp];
      
  }

    

    _S8_lin = (2.0 * k_2fs * _c8[qp] * (_I8[qp]) +
                  4.0 * k_1fs/(2 * k_2fs) * std::exp(k_2fs * std::pow(_I8[qp], 2.0))) * _fXs[qp].contract(_CLin) * _fXs[qp] ;
    



     
   
    _S_lin =  _S1_lin + _S4_lin + _S6_lin + _S8_lin;

   stressLin = _F[qp] * _S_lin + H * _S[qp];

}

void
CompressibleHolzapfel::compute_C_lin(unsigned const &qp, RealTensorValue const &H)
{
   _CLin = (_F[qp].transpose() * H + H.transpose() * _F[qp]);
}

 
