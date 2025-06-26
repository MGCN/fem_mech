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

#include "GuccioneCostaRot.h"

registerMooseObject("MECHApp", GuccioneCostaRot);

template <>
InputParameters
validParams<GuccioneCostaRot>()
{
  InputParameters params = validParams<ElasticMaterial>();
  params.addRequiredParam<Real>("mu", "mu");
  params.addRequiredParam<Real>("bf", "bf");
  params.addRequiredParam<Real>("bt", "bt");
  params.addRequiredParam<Real>("bfs", "bfs");
  // params.addParam<std::string>("base_name", "Material property base name");

  return params;
}

GuccioneCostaRot::GuccioneCostaRot(const InputParameters &parameters)
    : ElasticMaterial(parameters), mu(getParam<Real>("mu")),
      bf(getParam<Real>("bf")), bfs(getParam<Real>("bfs")),
      bt(getParam<Real>("bt")),
      _Ebar(declareProperty<RealTensorValue>("strain_bar")),
      _EbarRot(declareProperty<RealTensorValue>("strain_bar_rot")),
      _C(declareProperty<RealTensorValue>("CauchyGreen")),
      _invC(declareProperty<RealTensorValue>("inverseCauchyGreen")),
      _EtimesQ(declareProperty<RealTensorValue>("EtimesQ")),
      _SbarRot(declareProperty<RealTensorValue>("SbarRot")),
      _Sbar(declareProperty<RealTensorValue>("Sbar")),
      _S(declareProperty<RealTensorValue>("S")),
      _Pdev(declareProperty<RealTensorValue>("Pdev")),
      _R(getMaterialProperty<RealTensorValue>("rotationTensor")),
      _J23(declareProperty<Real>("J23")),
      _stiffening(declareProperty<Real>("stiffening"))
{
  _Q = RealTensorValue(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
  _Q(0, 0) = bf;
  _Q(0, 1) = bfs;
  _Q(0, 2) = bfs;
  _Q(1, 0) = bfs;
  _Q(2, 0) = bfs;
  _Q(1, 1) = bt;
  _Q(1, 2) = bt;
  _Q(2, 1) = bt;
  _Q(2, 2) = bt;

 _Q(0, 0) = 18.5;
  _Q(0, 1) = 2.8;
  _Q(0, 2) = 2.8;
  _Q(1, 0) = 2.8;
  _Q(2, 0) = 2.8;
  _Q(1, 1) = 3.58;
  _Q(1, 2) = 2.8;
  _Q(2, 1) = 2.8;
  _Q(2, 2) = 3.58;


}

void
GuccioneCostaRot::computeQpPropertiesDerived()
{
  // here you have available _F, _U, _stress, J because GuccioneCosta derives
  // from ElasticMaterial

  RealTensorValue const &F = _F[_qp];

  _C[_qp] = F.transpose() * F;

  RealTensorValue const &C = _C[_qp];

  _invC[_qp] = C.inverse();

  RealTensorValue invC = _invC[_qp];

  _J23[_qp] = std::pow(_J[_qp], -2. / 3.);

  Real const &J23 = _J23[_qp];

  _Ebar[_qp] = 0.5 * (J23 * C - _Id);

  _EbarRot[_qp] = _R[_qp].transpose() * _Ebar[_qp] * _R[_qp];

  RealTensorValue const &EbarRot = _EbarRot[_qp];

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      (_EtimesQ[_qp])(i, j) = _Q(i, j) * EbarRot(i, j);

  RealTensorValue const &EtimesQ = _EtimesQ[_qp];

  _stiffening[_qp] = mu * std::exp(EtimesQ.contract(EbarRot));

  Real const &st = _stiffening[_qp];

  _SbarRot[_qp] = st * EtimesQ;

  _Sbar[_qp] = _R[_qp] * _SbarRot[_qp] * _R[_qp].transpose();

  RealTensorValue const &Sbar = _Sbar[_qp];

  _S[_qp] = J23 * (Sbar - (1.0 / 3.0) * Sbar.contract(C) * invC);

  _P[_qp] += F * _S[_qp];
    
  _Pdev[_qp] = _P[_qp];
}

void
GuccioneCostaRot::evaluate_stress_lin(unsigned const &qp,
                                      RealTensorValue const &H,
                                      RealTensorValue &stressLin)
{
  compute_Ebar_prime(qp, H);
  _EbarH_Rot = _R[qp].transpose() * _EbarH * _R[qp];

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      _EH_Times_Q(i, j) = _EbarH_Rot(i, j) * _Q(i, j);

  Real const &q = _stiffening[qp];

  RealTensorValue temp =
      2.0 * _EtimesQ[qp] * (_EtimesQ[qp].contract(_EbarH_Rot)) + _EH_Times_Q;

  _SbarH = q * _R[qp] * temp * _R[qp].transpose();

  stressLin = -1.0 / 3.0 * _invC[qp].contract(_CH) * _Pdev[qp] + H * _S[qp] +
              _J23[qp] * _F[qp] *
                  (_SbarH - 1.0 / 3.0 * _SbarH.contract(_C[qp]) * _invC[qp] -
                   1.0 / 3.0 * _Sbar[qp].contract(_CH) * _invC[qp] +
                   1.0 / 3.0 * _Sbar[qp].contract(_C[qp]) * _invC[qp] * _CH *
                       _invC[qp]);
}

void
GuccioneCostaRot::compute_Ebar_prime(unsigned const &qp,
                                     RealTensorValue const &H)
{

  _CH = _F[qp].transpose() * H + H.transpose() * _F[qp];
  _EbarH = 0.5 * _J23[qp] * (_CH - 1.0 / 3.0 * (_CH * _invC[qp]).tr() * _C[qp]);
}
