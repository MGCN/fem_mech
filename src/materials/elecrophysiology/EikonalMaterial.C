#include "EikonalMaterial.h"
//#include "TensorHelpers.h"

registerMooseObject("MECHApp", EikonalMaterial);

template <>
InputParameters
validParams<EikonalMaterial>()
{
  InputParameters params = validParams<Material>();
    params.addRequiredParam<RealVectorValue>("sigma_i","forgot intracellular diffusion coefficients");
    //params.addRequiredParam<RealVectorValue>("sigma_e","forgot extracellular diffusion coefficients");
    //params.addRequiredParam<Real>("C_m","forgot membrane conductance");
    //params.addRequiredParam<Real>("Chi","forgot surface per volume");

  return params;
}

EikonalMaterial::EikonalMaterial(
    const InputParameters &parameters)
    : Material(parameters),
      _fXf(getMaterialProperty<RealTensorValue>("fiberOuterFiber")),
      _sXs(getMaterialProperty<RealTensorValue>("sheetOuterSheet")),
      _nXn(getMaterialProperty<RealTensorValue>("normalOuterNormal")),
      sigma_i(getParam<RealVectorValue>("sigma_i")),
      //sigma_e(getParam<RealVectorValue>("sigma_e")),
      //_mono(getParam<RealVectorValue>("mono_conductivity")),
      //C_m(getParam<Real>("C_m")),
      //Chi(getParam<Real>("Chi")),
      _Kmono(declareProperty<RealTensorValue>("_K"))
      //time_coefficient(declareProperty<Real>("time_coefficient"))
      
{}

void
EikonalMaterial::computeQpProperties()
{
    //RealVectorValue _mono;
    
  /*for (unsigned int i=0; i<3; ++i)
  {
    if (sigma_i(i)<1e-10 && sigma_e(i)<1e-10)
      mooseError("Both conductivies are zeros");
      
    _mono(i)=(sigma_i(i) * sigma_e(i))/(sigma_i(i) + sigma_e(i)) ;
  }*/

  _Kmono[_qp] = sigma_i(0) * _fXf[_qp] +
                sigma_i(1) * _sXs[_qp] +
                sigma_i(2) * _nXn[_qp];


    //time_coefficient[_qp] = C_m * Chi;

}


