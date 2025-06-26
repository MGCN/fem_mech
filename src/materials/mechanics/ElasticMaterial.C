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

#include "ElasticMaterial.h"
#include "MooseMesh.h"

#include "libmesh/quadrature.h"

template <>
InputParameters
validParams<ElasticMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredCoupledVar("disp_x","disp_x");
  params.addRequiredCoupledVar("disp_y","disp_y");
  params.addCoupledVar("disp_z","disp_z");
  params.addCoupledVar("pressure","pressure");
    
  params.addParam<bool>("store_variable_older", false, "Parameter which indicates whether the older stress state, required for HHT time integration, needs to be stored");

  return params;
}

ElasticMaterial::ElasticMaterial(const InputParameters &parameters)
    : Material(parameters),
      _dim(_mesh.dimension()),
      _J(declareProperty<Real>("deformationDeterminant")),
      _Ce(declareProperty<RealTensorValue>("CauchyGreenElastic")),
      _Se(declareProperty<RealTensorValue>("SecondPiolaElastic")),
      _pressure_material(declareProperty<Real>("pressureFromMaterial")),
      _grad_disp_x(coupledGradient("disp_x")),
      _grad_disp_y(coupledGradient("disp_y")),
      _grad_disp_z(_mesh.dimension() == 3 ? coupledGradient("disp_z") : _grad_zero),
      _has_pressure(parameters.isParamValid("pressure") ? 1 : 0),
      _pressure(parameters.isParamValid("pressure") ? coupledValue("pressure"):_zero),
      _U(declareProperty<RealTensorValue>("displacementGradient")),
      _F(declareProperty<RealTensorValue>("deformationGradient")),
      _invFtr(declareProperty<RealTensorValue>("invFtr")),
      _P(declareProperty<RealTensorValue>("firstPiolaKirchhofTensor")),
      _stress(declareProperty<RealTensorValue>("CauchyStress")),
      materialPointer(declareProperty<ElasticMaterial *>("materialPointer"))
{
//   if (_store_variable_older)
//   {
//      declarePropertyOlder<RealTensorValue>("invFtr");
//      declarePropertyOld<RealTensorValue>("invFtr");
//      declarePropertyOld<Real>("deformationDeterminant");
//      declarePropertyOlder<Real>("deformationDeterminant");
//      declarePropertyOld<RealTensorValue>("firstPiolaKirchhofTensor");
//      declarePropertyOlder<RealTensorValue>("firstPiolaKirchhofTensor");
//   }

  _Id = RealTensorValue(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);

  for (unsigned i = 0; i < _dim; ++i)
  {
    _Id(i, i) = 1.0;
  }
}

void
ElasticMaterial::initQpStatefulProperties()
{
  //_J[_qp] = 1.0;
  //_P[_qp]=RealTensorValue(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
  
}

void
ElasticMaterial::computeQpProperties()
{
    
  materialPointer[_qp] = this;
    
  _pressure_material[_qp]=_pressure[_qp];
  

    if (_dim == 2)
    {
        _U[_qp]=RealTensorValue(_grad_disp_x[_qp](0),_grad_disp_x[_qp](1),0.0,
                                _grad_disp_y[_qp](0),_grad_disp_y[_qp](1),0.0,
                                0.0                 ,0.0                 ,0.0);
        
        _F[_qp]=_U[_qp]+_Id;
        
        _J[_qp] = _F[_qp](0, 0) * _F[_qp](1, 1) - _F[_qp](1, 0) * _F[_qp](0, 1);
        
        RealTensorValue Ftemp = _F[_qp];

        Ftemp(2, 2) = 1.0;
        
        _invFtr[_qp] = Ftemp.transpose().inverse();
       
        for (int i = 0; i < 3; ++i)
        {
            _invFtr[_qp](i, 2) = 0.0;
            _invFtr[_qp](2, i) = 0.0;
        }
    }
    else if (_dim == 3)
    {
        _U[_qp] = RealTensorValue(_grad_disp_x[_qp], _grad_disp_y[_qp], _grad_disp_z[_qp]);
        _F[_qp] = _U[_qp] + _Id;
        _J[_qp] = _F[_qp].det();
        _invFtr[_qp] = _F[_qp].transpose().inverse();

   }

    // This can be removed when all materials will be implemented with _P= instead of _P+=
    _P[_qp] = RealTensorValue(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    
    computeQpPropertiesDerived();
    
//    for (int ii=0; ii<3; ++ii)
//        for (int jj=0; jj<3; ++jj)
//        {
//            _Pdev[_qp](ii,jj)=_P[_qp](ii,jj);
//        }
//    
//    if (_has_pressure)
//        for (int ii=0; ii<3; ++ii)
//            for (int jj=0; jj<3; ++jj)
//            {
//                _P[_qp](ii,jj)=_Pdev[_qp](ii,jj)+_pressure[_qp]*_J[_qp]*_invFtr[_qp](ii,jj);
//            }
    
    if(_has_pressure){
        _P[_qp]+=_pressure[_qp]*_J[_qp]*_invFtr[_qp];
      }
    
   _stress[_qp]=_P[_qp]*_F[_qp].transpose()/_J[_qp];
    
    
}
