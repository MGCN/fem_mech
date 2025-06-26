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
#include "NeoHookean.h"
registerMooseObject("MECHApp", NeoHookean);

template <>

InputParameters
validParams<NeoHookean>()
{
  InputParameters params = validParams<ElasticMaterial>();
  params.addRequiredParam<Real>("mu", "mu");
  params.addRequiredParam<Real>("lambda", "lambda");
  return params;
}

NeoHookean::NeoHookean(const InputParameters &parameters)
    : ElasticMaterial(parameters), mu(getParam<Real>("mu")),
      lambda(getParam<Real>("lambda")),
      _C(declareProperty<RealTensorValue>("CauchyGreen")),
      _invC(declareProperty<RealTensorValue>("inverseCauchyGreen"))
{
}

void
NeoHookean::computeQpPropertiesDerived()
{
  // here you have available _F, _U, _stress, J because NeoHookean derives
  // from ElasticMaterial

    RealTensorValue const &F = _F[_qp];

    Real const &J = _J[_qp];


     
    
    if (_mesh.dimension() == 3)
    {
        _C[_qp]      = F.transpose() * F;
        
        _invC[_qp]   = _C[_qp].inverse();

    }



    _P[_qp] = mu * (F - _invFtr[_qp]) + lambda * std::log(_J[_qp]) * _invFtr[_qp];
}

void
NeoHookean::evaluate_stress_lin(unsigned const &qp, RealTensorValue const &H,
                                RealTensorValue &stressLin)

{

   Real const &J = _J[qp];

   RealTensorValue const &invFtr = _invFtr[qp];

   RealTensorValue const &F = _F[_qp];

   RealTensorValue const &invCe = _invC[qp];

   RealTensorValue const &Ce = _C[_qp];

   RealTensorValue  CLin = 1.0 * (F.transpose() * H + H.transpose() * F);

   RealTensorValue  S = mu * _Id + (lambda * std::log(_J[_qp]) - 1.0 * mu ) * invCe;

   //RealTensorValue  Slin = - 1.0 * (1.0 * lambda * std::log(J) - 1.0 * mu) * invCe * CLin * invCe + lambda * 1.0/2.0 * (invCe.contract(CLin)) * invCe;

   //RealTensorValue  Plin = H * S + Slin * F;
   
   RealTensorValue Plin = mu * H - (1.0 * lambda * std::log(J) - 1.0 * mu) * invFtr * H.transpose() * invFtr + lambda * invFtr.contract(H) * invFtr;

  stressLin = Plin;
}
