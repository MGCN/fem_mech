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
/****************************************************************/


#include "ActiveStressFunction_timeactivated.h"

#include "Material.h"
#include "Function.h"

registerMooseObject("MECHApp", ActiveStressFunction_timeactivated);

template <>
InputParameters
validParams<ActiveStressFunction_timeactivated>()
{

  // inherit the parameters of the Kernels:
  InputParameters params = validParams<Kernel>();

  // specify for which component to solve:
  params.addRequiredParam<unsigned>("component", "component");
    
  //params.addRequiredParam<unsigned>("function_mod", "function_mod");
    params.addRequiredParam<Real>("function_mod", "function_mod");

//time activation
    params.addRequiredCoupledVar(
        "activationTime",
        "Aux variable for time activation given by Eikonal model");

 //specify if time map is in second (true if second false if millisecond)
 params.addRequiredParam<bool>("time_unit", "time_unit");
 
  return params;
}

ActiveStressFunction_timeactivated::ActiveStressFunction_timeactivated(const InputParameters &parameters)
    :

      // inherit the parameters of the Kernels:
      Kernel(parameters),

      // specify for which component to solve:
      _component(getParam<unsigned>("component")),

      // inherit some material properties:
      _fXf(getMaterialProperty<RealTensorValue>("fiberOuterFiber")),
      _F(getMaterialProperty<RealTensorValue>("deformationGradient")),

      // the driving function
      //_activeFunction_modulus(getFunction("function_1")),
     _activeFunction_modulus(getParam<Real>("function_mod")),

     // time activation map
     _time_activation_map(coupledValue("activationTime")),

     //bool time_unit
     _time_unit(getParam<bool>("time_unit"))
{
}

Real
ActiveStressFunction_timeactivated::computeQpResidual()
{

  // generate the test function
  V = RealTensorValue(zeros, zeros, zeros);

  for (unsigned k = 0; k < 3; ++k)
    V(_component, k) = _grad_test[_i][_qp](k);

  RealTensorValue Pa;

  // evaluation of the function at the current time for the quadrature point
Real factor;
  /*Real factor = 0.0*(_t)*(_t<=0.2 + _time_activation_map.value(_t,_q_point[_qp]))+_activeFunction_modulus.value(_t, _q_point[_qp])*(_t - (0.2 + ( _time_activation_map.value(_t, _q_point[_qp])/1000)))*(_t>0.2 + (_time_activation_map.value(_t, _q_point[_qp])/1000));*/
    if (_time_unit){
	/*factor = 0.0*(_t)*(_t<=0.0001 + _time_activation_map.value(_t, _q_point[_qp]))+_activeFunction_modulus*(_t - (0.0001 + ( _time_activation_map.value(_t, _q_point[_qp]))))*(_t>0.0001 + (_time_activation_map.value(_t, _q_point[_qp])));*/
        
	    /*factor = 0.0*(_t)*(_t<=_time_activation_map[_qp]/1.0)+(_activeFunction_modulus/(0.2*0.2))*(_t -(_time_activation_map[_qp]/1.0))*(_t -(_time_activation_map[_qp]/1.0))*(_t>(_time_activation_map[_qp]/1.0));*/

	            factor = 0.0*(_t)*(_t<=_time_activation_map[_qp])+(_activeFunction_modulus/((0.2-_time_activation_map[_qp])*(0.2- _time_activation_map[_qp])*(0.2- _time_activation_map[_qp])))*(_t -_time_activation_map[_qp])*(_t -_time_activation_map[_qp])*(_t -_time_activation_map[_qp])*(_t>(_time_activation_map[_qp]));

       /* if (factor != 0){
            std::cout<<"factor:"<<factor<<std::endl;*/
        }

    
    else {
	   /*factor = 0.0*(_t)*(_t<=0.0001 + _time_activation_map.value(_t,_q_point[_qp])/1000)+_activeFunction_modulus*(_t - (0.0001 + ( _time_activation_map.value(_t, _q_point[_qp])/1000)))*(_t>0.0001 + (_time_activation_map.value(_t, _q_point[_qp])/1000));*/
        factor = 0.0*(_t)*(_t<=_time_activation_map[_qp]/1000)+_activeFunction_modulus*(_t -(_time_activation_map[_qp]/1000))*(_t>(_time_activation_map[_qp]/1000));
    }
  RealTensorValue F = _F[_qp];

  // computation of the active tension
  Pa = factor * F * _fXf[_qp];

  // contraction of the active tension with the test function
  return Pa.contract(V);
}

Real
ActiveStressFunction_timeactivated::computeQpJacobian()
{

  // Computation of the diagonal entries of the Jacobian matrix

  H = RealTensorValue(zeros, zeros, zeros);
  V = RealTensorValue(zeros, zeros, zeros);

  for (unsigned k = 0; k < 3; ++k)
  {
    H(_component, k) = _grad_phi[_j][_qp](k);
  }

  for (unsigned k = 0; k < 3; ++k)
    V(_component, k) = _grad_test[_i][_qp](k);

  Real factor;
  // evaluation of the function at the current time for the quadrature point
    if (_time_unit){
    /*factor = 0.0*(_t)*(_t<=0.0001 + _time_activation_map.value(_t, _q_point[_qp]))+_activeFunction_modulus*(_t - (0.0001 + ( _time_activation_map.value(_t, _q_point[_qp]))))*(_t>0.0001 + (_time_activation_map.value(_t, _q_point[_qp])));*/
        
	    /*factor = 0.0*(_t)*(_t<=_time_activation_map[_qp])+(_activeFunction_modulus/(0.2*0.2))*(_t -(_time_activation_map[_qp]/1.0))*(_t -(_time_activation_map[_qp]/1.0))*(_t>(_time_activation_map[_qp]));*/

	                      factor = 0.0*(_t)*(_t<=_time_activation_map[_qp])+(_activeFunction_modulus/((0.2-_time_activation_map[_qp])*(0.2- _time_activation_map[_qp])*(0.2- _time_activation_map[_qp])))*(_t -_time_activation_map[_qp])*(_t -_time_activation_map[_qp])*(_t -_time_activation_map[_qp])*(_t>(_time_activation_map[_qp]));
						      }

    else {
       /*factor = 0.0*(_t)*(_t<=0.0001 + _time_activation_map.value(_t,_q_point[_qp])/1000)+_activeFunction_modulus*(_t - (0.0001 + ( _time_activation_map.value(_t, _q_point[_qp])/1000)))*(_t>0.0001 + (_time_activation_map.value(_t, _q_point[_qp])/1000));*/
        factor = 0.0*(_t)*(_t<=_time_activation_map[_qp]/1000)+_activeFunction_modulus*(_t -(_time_activation_map[_qp]/1000))*(_t>(_time_activation_map[_qp]/1000));
    }

  // computation of the linearization of the active tension
  RealTensorValue Pa = factor * H * _fXf[_qp];

  // contraction of the linearization of the active tension with the test
  // function
  return Pa.contract(V);
}

Real
ActiveStressFunction_timeactivated::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Computation of the off diagonal entries of the Jacobian matrix 
  // jvar is the variable w.r.t. which we are derivating (the considered coloumn). 

  if (jvar > 2)
    return 0.0;

  H = RealTensorValue(zeros, zeros, zeros);
  V = RealTensorValue(zeros, zeros, zeros);

  for (unsigned k = 0; k < 3; ++k)
  {
    H(jvar, k) = _grad_phi[_j][_qp](k);
  }

  for (unsigned k = 0; k < 3; ++k)
    V(_component, k) = _grad_test[_i][_qp](k);
Real factor;
    if (_time_unit){
    /*factor = 0.0*(_t)*(_t<=0.0001 + _time_activation_map.value(_t, _q_point[_qp]))+_activeFunction_modulus*(_t - (0.0001 + ( _time_activation_map.value(_t, _q_point[_qp]))))*(_t>0.0001 + (_time_activation_map.value(_t, _q_point[_qp])));*/
       /* factor = 0.0*(_t)*(_t<=_time_activation_map[_qp])+(_activeFunction_modulus/(0.2*0.2))*(_t -(_time_activation_map[_qp]/1.0))*(_t -(_time_activation_map[_qp]/1.0))*(_t>(_time_activation_map[_qp]));*/
                            factor = 0.0*(_t)*(_t<=_time_activation_map[_qp])+(_activeFunction_modulus/((0.2-_time_activation_map[_qp])*(0.2- _time_activation_map[_qp])*(0.2- _time_activation_map[_qp])))*(_t -_time_activation_map[_qp])*(_t -_time_activation_map[_qp])*(_t -_time_activation_map[_qp])*(_t>(_time_activation_map[_qp]));
						    }
    
    else {
       /*factor = 0.0*(_t)*(_t<=0.0001 + _time_activation_map.value(_t,_q_point[_qp])/1000)+_activeFunction_modulus*(_t - (0.0001 + ( _time_activation_map.value(_t, _q_point[_qp])/1000)))*(_t>0.0001 + (_time_activation_map.value(_t, _q_point[_qp])/1000));*/
        factor = 0.0*(_t)*(_t<=_time_activation_map[_qp]/1000)+_activeFunction_modulus*(_t -(_time_activation_map[_qp]/1000))*(_t>(_time_activation_map[_qp]/1000));
    }
  // linearization of the active tension
  RealTensorValue Pa = factor * H * _fXf[_qp];

  // contraction of the linearized active tension with the test function
  return Pa.contract(V);
}
