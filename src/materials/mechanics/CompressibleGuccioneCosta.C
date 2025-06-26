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
/*                Prepared by Marco Favino,                     */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/* Construction of Guccione Costa Material.                     */
/****************************************************************/


#include "CompressibleGuccioneCosta.h"

registerMooseObject("MECHApp", CompressibleGuccioneCosta);

template <>
InputParameters
validParams<CompressibleGuccioneCosta>()
{
  InputParameters params = validParams<ElasticMaterial>();
  params.addRequiredParam<Real>("mu", "mu");
  params.addRequiredParam<Real>("bf", "bf");
  params.addRequiredParam<Real>("bt", "bt");
  params.addRequiredParam<Real>("bfs", "bfs");

  return params;
}

CompressibleGuccioneCosta::CompressibleGuccioneCosta(
    InputParameters const &parameters)
    : ElasticMaterial(parameters), mu(getParam<Real>("mu")),
      bf(getParam<Real>("bf")), bfs(getParam<Real>("bfs")),
      bt(getParam<Real>("bt")), _E(declareProperty<RealTensorValue>("strain")),
      _C(declareProperty<RealTensorValue>("CauchyGreen")),
      _EtimesQ(declareProperty<RealTensorValue>("EtimesQ")),
      _S(declareProperty<RealTensorValue>("S")),
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
}

void
CompressibleGuccioneCosta::computeQpPropertiesDerived()
{
  // here you have available _F, _U, _stress, J because GuccioneCosta derives
  // from ElasticMaterial

  RealTensorValue const &F = _F[_qp];

  _C[_qp] = F.transpose() * F;

  RealTensorValue const &C = _C[_qp];

  _E[_qp] = 0.5 * (C - _Id);

  RealTensorValue const &E = _E[_qp];

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      (_EtimesQ[_qp])(i, j) = _Q(i, j) * E(i, j);

  RealTensorValue const &EtimesQ = _EtimesQ[_qp];

  _stiffening[_qp] = mu * std::exp(EtimesQ.contract(E));

  _S[_qp] = _stiffening[_qp] * EtimesQ;

  RealTensorValue const &S = _S[_qp];

  _P[_qp] = F * S;

}

void
CompressibleGuccioneCosta::evaluate_stress_lin(unsigned const &qp,
                                               RealTensorValue const &H,
                                               RealTensorValue &stressLin)
{
  compute_E_prime(qp, H);

  RealTensorValue SH;

  RealTensorValue EHROT;

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      EHROT(i, j) = _EH(i, j) * _Q(i, j);

  Real const &q = _stiffening[qp];

  SH = q * (EHROT + 2.0 * _EtimesQ[qp].contract(_EH) * _EtimesQ[qp]);

  stressLin = _F[qp] * SH + H * _S[qp];
}

void
CompressibleGuccioneCosta::compute_E_prime(unsigned const &qp,
                                           RealTensorValue const &H)
{
  _EH = 0.5 * (_F[qp].transpose() * H + H.transpose() * _F[qp]);
}
