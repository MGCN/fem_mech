#include "MonodomainConductivity.h"
//#include "TensorHelpers.h"

registerMooseObject("MECHApp", MonodomainConductivity);

template <>
InputParameters
validParams<MonodomainConductivity>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<std::vector<Real>>(
      "extraCellularConductivities",
      "Conductivities in fibre, sheet and sheet-normal direction. Unit: S/cm.");
    params.addRequiredParam<std::vector<Real>>(
      "intraCellularConductivities",
      "Conductivities in fibre, sheet and sheet-normal direction. Unit: S/cm.");

  return params;
}

MonodomainConductivity::MonodomainConductivity(
    const InputParameters &parameters)
    : Material(parameters),
      _fXf(getMaterialProperty<RealTensorValue>("fiberOuterFiber")),
      _sXs(getMaterialProperty<RealTensorValue>("sheetOuterSheet")),
      _nXn(getMaterialProperty<RealTensorValue>("normalOuterNormal")),
      _intraConductivities(getParam<std::vector<Real>>("intraCellularConductivities")),
      _extraConductivities(getParam<std::vector<Real>>("extraCellularConductivities")),
      _Gmono(declareProperty<RealTensorValue>("conductivity")),
      _Gintra(declareProperty<RealTensorValue>("Gintra"))
{
  if (_intraConductivities.size() != 3)
    mooseError("MonodomainConductivity: conductivities must contain exactly 3 "
               "numbers");
  if (_extraConductivities.size() != 3)
    mooseError("MonodomainConductivity: conductivities must contain exactly 3 "
               "numbers");

  _monoConductivities.resize(3);

  for (unsigned int i=0; i<3; ++i)
  {
    if (_intraConductivities[i]<1e-10 && _extraConductivities[i]<1e-10)
      mooseError("Both conductivies are zeros");
    _monoConductivities[i]=(_intraConductivities[i] * _extraConductivities[i])/ (_intraConductivities[i] + _extraConductivities[i]) ;
  }

}

void
MonodomainConductivity::computeQpProperties()
{


  _Gmono[_qp] = _monoConductivities[0] * _fXf[_qp] +
                _monoConductivities[1] * _sXs[_qp] +
                _monoConductivities[2] * _nXn[_qp];

  _Gintra[_qp] = _intraConductivities[0] * _fXf[_qp] +
                 _intraConductivities[1] * _sXs[_qp] +
                 _intraConductivities[2] * _nXn[_qp];


}

