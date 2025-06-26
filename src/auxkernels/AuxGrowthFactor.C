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
/*                Prepared by Maria Nestola,                    */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/*                                                              */
/*            Auxkernel to compute velocity field               */
/****************************************************************/

#include "AuxGrowthFactor.h"

registerMooseObject("MECHApp", AuxGrowthFactor);

template<>
InputParameters 
validParams<AuxGrowthFactor>()
{
    InputParameters params = validParams<AuxKernel>();
	params.addRequiredParam<Real>("kappa_plus", "kappa_plus");
	params.addRequiredParam<Real>("kappa_minus", "kappa_minus");
	params.addRequiredParam<Real>("theta_plus", "theta_plus");
	params.addRequiredParam<Real>("theta_minus", "theta_minus");   
	params.addRequiredParam<Real>("m_plus", "m_plus");
	params.addRequiredParam<Real>("m_minus", "m_minus");  
	params.addRequiredCoupledVar("threshold","threshold");
    return params;
}
AuxGrowthFactor::AuxGrowthFactor(const InputParameters & parameters) :
AuxKernel(parameters),
// _u_old(valueOld()),
_k_plus(getParam<Real>("kappa_plus")),
_k_minus(getParam<Real>("kappa_minus")),
_threshold(coupledValue("threshold")),
_theta_plus(getParam<Real>("theta_plus")),
_theta_minus(getParam<Real>("theta_minus")),
_m_plus(getParam<Real>("m_plus")),
_m_minus(getParam<Real>("m_minus")),
u_old(uOld())

{
}

Real
AuxGrowthFactor::computeValue()
{
    
      Real k_theta = 0;
      
      if(_threshold[_qp] > -1.0e-12)
      {
        
        Real g_f = _k_plus * ( _theta_plus  - u_old[_qp]) / (_theta_plus - 1.0);

        k_theta = std::pow(g_f, _m_plus) * _threshold[_qp];
      }
      else{
        
        Real g_f = _k_minus * ( u_old[_qp]  - _theta_minus) / (1.0 - _theta_minus);

        k_theta =  std::pow(g_f, _m_minus) * _threshold[_qp];
      }

      
      if(_t_step==1) k_theta=1.0;

      std::cout<<"g_f_minus"<<k_theta<<std::endl;

      return k_theta;
//vel_old + ( _dt * ( 1 -_gamma)) * _accel_old[_qp] + _gamma * _dt * _accel[_qp];
}
