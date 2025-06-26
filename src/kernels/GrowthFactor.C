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
/* Kernel to ensure the mass balance. To be used only coupled   */
/* with IncompressibleStressDivergenceLibmeshTensor             */
/****************************************************************/

#include "GrowthFactor.h"
#include "MooseMesh.h"

#include "libmesh/quadrature.h"

registerMooseObject("MECHApp", GrowthFactor);

template <>
InputParameters
validParams<GrowthFactor>()
{

  // inherit the parameters of the Kernels:
  InputParameters params = validParams<Kernel>();


  params.addRequiredParam<Real>("kappa_plus", "kappa_plus");
  params.addRequiredParam<Real>("kappa_minus", "kappa_minus");
  params.addRequiredParam<Real>("theta_plus", "theta_plus");
  params.addRequiredParam<Real>("theta_minus", "theta_minus");   
  params.addRequiredParam<Real>("m_plus", "m_plus");
  params.addRequiredParam<Real>("m_minus", "m_minus");  
  params.addRequiredCoupledVar("threshold","threshold");

  return params;
}

GrowthFactor::GrowthFactor(InputParameters const &parameters)
    :

      // inherit the parameters of the Kernels:
      Kernel(parameters),
      _k_plus(getParam<Real>("kappa_plus")),
      _k_minus(getParam<Real>("kappa_minus")),
      _threshold(coupledValue("threshold")),
      _theta_plus(getParam<Real>("theta_plus")),
      _theta_minus(getParam<Real>("theta_minus")),
      _m_plus(getParam<Real>("m_plus")),
      _m_minus(getParam<Real>("m_minus"))

{

}

Real
GrowthFactor::computeQpResidual()
{
     
      Real k_theta = 0;
      
      if(_threshold[_qp] > -1.0e-12)
      {
        Real g_f = _k_plus * ( _theta_plus  - _u[_qp]) / (_theta_plus - 1.0);

        //std::cout<<"g_f_plus"<<g_f<<std::endl;

        k_theta = std::pow(g_f, _m_plus) * _threshold[_qp];
      }
      else{
        
        Real g_f = _k_minus * ( _u[_qp]  - _theta_minus) / (1.0 - _theta_minus);

        //std::cout<<"g_f_minus"<<g_f<<std::endl;

        k_theta =  std::pow(g_f, _m_minus) * _threshold[_qp];
      }

      return k_theta * _test[_i][_qp];
}

Real
GrowthFactor::computeQpOffDiagJacobian(unsigned int jvar)
{
    return 0.0;
}

Real
GrowthFactor::computeQpJacobian()
{

      Real k_theta     = 0;

      Real k_theta_lin = 0;

      if(_threshold[_qp] > -1.0e-12)
      {
        Real g_f = _k_plus * ( _theta_plus  - _u[_qp]) / (_theta_plus - 1.0);

        k_theta = std::pow(g_f, _m_plus) * _threshold[_qp];

        k_theta_lin =  k_theta * _m_plus/(_u[_qp] - _theta_plus);
      }
      else{
        
        Real g_f = _k_minus * ( _u[_qp]  - _theta_minus) / (1.0 - _theta_minus);

        k_theta =  std::pow(g_f, _m_minus) * _threshold[_qp];

        k_theta_lin =  k_theta * _m_minus/(_u[_qp] - _theta_minus);

      }

      return k_theta_lin * _phi[_j][_qp] * _test[_i][_qp];
}
