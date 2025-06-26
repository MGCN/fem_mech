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
/* Construction of Guccione Costa Material.                     */
/****************************************************************/

#include "GuccioneCosta.h"
#include "libmesh/quadrature.h"

#include "MooseMesh.h"


registerMooseObject("MECHApp", GuccioneCosta);


template <>
InputParameters
validParams<GuccioneCosta>()
{
  InputParameters params = validParams<ElasticMaterial>();
  params.addRequiredParam<Real>("mu", "mu");
  params.addRequiredParam<Real>("bf", "bf");
  params.addRequiredParam<Real>("bt", "bt");
  params.addRequiredParam<Real>("bfs", "bfs");
  return params;
}

GuccioneCosta::GuccioneCosta(const InputParameters &parameters)
    : ElasticMaterial(parameters), mu(getParam<Real>("mu")),
      bf(getParam<Real>("bf")), bfs(getParam<Real>("bfs")),
      bt(getParam<Real>("bt")),
      _Ebar(declareProperty<RealTensorValue>("incompressibleStrain")),
      _C(declareProperty<RealTensorValue>("CauchyGreen")),
      _invC(declareProperty<RealTensorValue>("inverseCauchyGreen")),
      _EbarTimesQ(declareProperty<RealTensorValue>("EbarTimesQ")),
      _Sbar(declareProperty<RealTensorValue>("Sbar")),
      _S(declareProperty<RealTensorValue>("S")),
      _Pdev(declareProperty<RealTensorValue>("Pdev")),
      _J23(declareProperty<Real>("J23")),
      _stiffening(declareProperty<Real>("stiffening"))
{
  _Q = RealTensorValue(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
  _Q(0, 0) = bf; //fibres
  _Q(0, 1) = bfs;
  _Q(0, 2) = bfs;
  _Q(1, 0) = bfs;
  _Q(2, 0) = bfs;
  _Q(1, 1) = bt; //sheet
  _Q(1, 2) = bt;
  _Q(2, 1) = bt;
  _Q(2, 2) = bt;
}

void
GuccioneCosta::computeQpPropertiesDerived()
{
  // here you have available _F, _U, _stress, J because GuccioneCosta derives
  // from ElasticMaterial


  RealTensorValue const &F = _F[_qp];

  _C[_qp] = F.transpose() * F;

  if(_dim==2)
    _C[_qp](2,2)=1.0;

  RealTensorValue const &C = _C[_qp];


  _invC[_qp] = C.inverse();

   if (_dim == 2) {

        const double & one= 1.0;
      
        RealTensorValue Ctemp = C;

        Ctemp(2, 2) = one;
        _invC[_qp] = Ctemp.inverse();
        
        for (int i = 0; i < 3; ++i)
        {
          _invC[_qp](i, 2) = 0.0;
        }
  }
   
  if (_dim == 3) {
    _invC[_qp] = C.inverse();
  
  }


  RealTensorValue invC = _invC[_qp];

  _J23[_qp] = std::pow(_J[_qp], -2. / _dim);
  Real const &J23 = _J23[_qp];



  _Ebar[_qp] = 0.5 * (J23 * C - _Id);
  RealTensorValue const &Ebar = _Ebar[_qp];

  for (int i = 0; i < _dim; ++i)
    for (int j = 0; j < _dim; ++j)
      (_EbarTimesQ[_qp])(i, j) = _Q(i, j) * Ebar(i, j);

  RealTensorValue const &EbarTimesQ = _EbarTimesQ[_qp];

  _stiffening[_qp] = mu * std::exp(EbarTimesQ.contract(Ebar));

  _Sbar[_qp] = _stiffening[_qp] * EbarTimesQ;
  RealTensorValue const &Sbar = _Sbar[_qp];
    
  _S[_qp] = J23 * (Sbar - (1.0 / _dim) * Sbar.contract(C) * invC);

  _P[_qp] += F * _S[_qp];
    
  _Pdev[_qp] = _P[_qp];
    
}

void
GuccioneCosta::compute_Ebar_prime(unsigned const &qp, RealTensorValue const &H)
{


  _CH = _F[qp].transpose() * H + H.transpose() * _F[qp];
  _EbarH = 0.5 * _J23[qp] * (_CH - 1.0 / _dim * (_CH * _invC[qp]).tr() * _C[qp]);
}

void
GuccioneCosta::evaluate_stress_lin(unsigned const &qp, RealTensorValue const &H,
                                   RealTensorValue &stressLin)
{


  compute_Ebar_prime(qp, H);

  for (int i = 0; i < _dim; ++i)
    for (int j = 0; j < _dim; ++j)
      (tempPrime)(i, j) = _Q(i, j) * _EbarH(i, j);

  _SbarH =
      _stiffening[qp] *
      (tempPrime + 2.0 * _EbarTimesQ[qp].contract(_EbarH) * _EbarTimesQ[qp]);

  stressLin = -1.0 / _dim * _invC[qp].contract(_CH) * _Pdev[qp] + H * _S[qp] +
              _J23[qp] * _F[qp] *
                  (_SbarH - 1.0 / _dim * _SbarH.contract(_C[qp]) * _invC[qp] -
                   1.0 / _dim * _Sbar[qp].contract(_CH) * _invC[qp] +
                   1.0 / _dim * _Sbar[qp].contract(_C[qp]) * _invC[qp] * _CH *
                       _invC[qp]);
}
