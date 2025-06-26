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

#include "LinearElastic.h"
#include "MooseMesh.h"

registerMooseObject("MECHApp", LinearElastic);

template <>
InputParameters
validParams<LinearElastic>()
{
  InputParameters params = validParams<ElasticMaterial>();
  params.addRequiredParam<Real>("mu", "mu");
  params.addRequiredParam<Real>("kappa", "kappa");
  return params;
}

LinearElastic::LinearElastic(InputParameters const &parameters)
    : ElasticMaterial(parameters),
      mu(getParam<Real>("mu")),
      kappa(getParam<Real>("kappa")),
      _trE(declareProperty<Real>("trE")),
      _E(declareProperty<RealTensorValue>("strain")),
      _S(declareProperty<RealTensorValue>("S"))

{
  if (kappa < 0.0)
    kappa = 0.0;

  // Computation of _lambda
  lambda = kappa - 2.0 * mu / _dim;
}

void
LinearElastic::computeQpPropertiesDerived()
{

    _F[_qp]=_Id;
    
    _invFtr[_qp]=_Id;

    _pressure_material[_qp]=0.0;
    
    RealTensorValue const &U = _U[_qp];
    
    _J[_qp]=1.0;
    
    _E[_qp] = 0.5 * (U.transpose() + U);
    RealTensorValue const &E = _E[_qp];

    _trE[_qp] = E.tr();
    Real const &trE = _trE[_qp];

    _stress[_qp] = 2.0 * mu * E + lambda * trE * _Id;

    // observe that here we have to replace
    _P[_qp] += _stress[_qp];
    
}

void
LinearElastic::evaluate_stress_lin(unsigned const &qp, RealTensorValue const &H,
                                   RealTensorValue &stressLin)
{
  RealTensorValue _strain_lin = 0.5 * (H + H.transpose());
  stressLin = 2.0 * mu * _strain_lin + lambda * _strain_lin.tr() * _Id;
}
