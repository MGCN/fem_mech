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
/*                  modified by Maria Nestola                     */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/* Fibers Material to compute fibers direction.                 */
/****************************************************************/

#include "FibersMaterial.h"

template <>
InputParameters
validParams<FibersMaterial>()
{
  InputParameters params = validParams<Material>();
  return params;
}

FibersMaterial::FibersMaterial(const InputParameters &parameters)
    : Material(parameters),
      _fiberDirection(declareProperty<RealVectorValue>("fiberDirection")),
      _gfiberDirection(declareProperty<RealVectorValue>("gfiberDirection")),
      _sheetDirection(declareProperty<RealVectorValue>("sheetDirection")),
      _normalDirection(declareProperty<RealVectorValue>("normalDirection")),
      _fXf(declareProperty<RealTensorValue>("fiberOuterFiber")),
      _gXg(declareProperty<RealTensorValue>("gfiberOutergFiber")),
      _sXs(declareProperty<RealTensorValue>("sheetOuterSheet")),
      _fXs(declareProperty<RealTensorValue>("fiberOuterSheet")),
      _nXn(declareProperty<RealTensorValue>("normalOuterNormal")),
      _rotationTensor(declareProperty<RealTensorValue>("rotationTensor"))
{
}

void
FibersMaterial::computeQpTensorProperties()
{
//  std::cout<<"_fiberDirection[_qp]"<<_fiberDirection[_qp]<<std::endl;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    {
      _fXf[_qp](i, j) = _fiberDirection[_qp](i)  * _fiberDirection[_qp](j);
      _gXg[_qp](i, j) = _gfiberDirection[_qp](i) * _gfiberDirection[_qp](j);
      _sXs[_qp](i, j) = _sheetDirection[_qp](i)  * _sheetDirection[_qp](j);
      _fXs[_qp](i, j) = _fiberDirection[_qp](i)  * _sheetDirection[_qp](j);
      _nXn[_qp](i, j) = _normalDirection[_qp](i) * _normalDirection[_qp](j);
    }

  // Here we symmetrize _fXs, maybe this operation is not needed or dangerous

  _fXs[_qp] = 0.5 * (_fXs[_qp] + _fXs[_qp].transpose());

// std::cout<<"f"<<_fiberDirection[_qp]<<std::endl;

// std::cout<<"s"<<_sheetDirection[_qp]<<std::endl;

// std::cout<<"n"<<_normalDirection[_qp]<<std::endl;

  for (int i = 0; i < 3; ++i)
  {
    _rotationTensor[_qp](i, 0) = _fiberDirection[_qp](i);
    _rotationTensor[_qp](i, 1) = _sheetDirection[_qp](i);
    _rotationTensor[_qp](i, 2) = _normalDirection[_qp](i);
  }

  // This part checks that _RotationTensor is orthonormal
  RealTensorValue OP = _rotationTensor[_qp].transpose() * _rotationTensor[_qp];

  // std::cout<<OP<<std::endl;

  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      if (i == j)
      {
        if (fabs(OP(i, j) - 1.0) > 1.0e-5)
        {
          //std::cout<<_rotationTensor[_qp](i, j)<<std::endl;
          mooseError("Rotation Tensor is not orthonormal.");
        }
      }
      else
      {
        if (fabs(OP(i, j)) > 1.0e-5)
        {
          //std::cout<<fabs(OP(i, j))<<std::endl;
          mooseError("Rotation Tensor is not orthonormal.");
        }

      }
    }
  }
}
