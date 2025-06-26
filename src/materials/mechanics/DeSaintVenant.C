
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


#include "DeSaintVenant.h"
#include "MooseMesh.h"

registerMooseObject("MECHApp", DeSaintVenant);

template <>
InputParameters
validParams<DeSaintVenant>()
{
  InputParameters params = validParams<ElasticMaterial>();
  params.addRequiredParam<Real>("mu", "mu");
  params.addRequiredParam<Real>("kappa", "kappa");
  params.addRequiredParam<Real>("epsilon", "epsilon");

  return params;
}

DeSaintVenant::DeSaintVenant(const InputParameters &parameters)
    : ElasticMaterial(parameters), mu(getParam<Real>("mu")),
      kappa(getParam<Real>("kappa")), epsilon(getParam<Real>("epsilon")),
      _E(declareProperty<RealTensorValue>("strain")),
      _C(declareProperty<RealTensorValue>("CauchyGreen")),
      _S(declareProperty<RealTensorValue>("S")),
      _Pdev(declareProperty<RealTensorValue>("Pdev"))
{

  if (kappa < epsilon)
    kappa = epsilon;

  // Computation of _lambda
  lambda = kappa - 2.0 * mu / _mesh.dimension();
}

void
DeSaintVenant::computeQpPropertiesDerived()
{
  // here you have available _F, _U, _stress, J because GuccioneCosta derives
  // from ElasticMaterial

  RealTensorValue const &F = _F[_qp];

  _C[_qp] = F.transpose() * F;

  RealTensorValue const &C = _C[_qp];

  _E[_qp] = 0.5 * (C - _Id);

  RealTensorValue const &E = _E[_qp];

  _S[_qp] = 2.0 * mu * E + lambda * E.tr() * _Id;

  _stress[_qp] = _S[_qp];

  _P[_qp] += F * _S[_qp];
  
  _Pdev[_qp]=_P[_qp];
}

void
DeSaintVenant::evaluate_stress_lin(unsigned const &qp, RealTensorValue const &H,
                                   RealTensorValue &stressLin)
{
  RealTensorValue _strain_lin =
      0.5 * (_F[qp].transpose() * H + H.transpose() * _F[qp]);
  RealTensorValue _S_lin =
      2.0 * mu * _strain_lin + lambda * _strain_lin.tr() * _Id;
  stressLin = _F[qp] * _S_lin + H * _S[qp];
}
