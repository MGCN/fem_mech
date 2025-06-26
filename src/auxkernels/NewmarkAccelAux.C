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
/*          Auxkernel to compute acceleration field             */
/****************************************************************/


#include "NewmarkAccelAux.h"

registerMooseObject("MECHApp", NewmarkAccelAux);


template<>
InputParameters validParams<NewmarkAccelAux>()
{
  InputParameters params = validParams<AuxKernel>();
    params.addRequiredCoupledVar("displacement","displacement variable");
    params.addRequiredCoupledVar("velocity","velocity variable");
    params.addRequiredParam<Real>("beta","beta parameter");
  return params;
}

NewmarkAccelAux::NewmarkAccelAux(const InputParameters & parameters) :
  AuxKernel(parameters),
   _disp_old(coupledValueOld("displacement")),
   _disp(coupledValue("displacement")),
   _vel_old(coupledValueOld("velocity")),
   _beta(getParam<Real>("beta")),
   u_old(uOld())

{
}

Real
NewmarkAccelAux::computeValue()
{
  if (!isNodal())
    mooseError("must run on a nodal variable");
    
  Real accel_old = u_old[_qp];
    
  if (_dt == 0)
    return accel_old;
    
  //Calculates acceeleration using Newmark time integration method
  return 1/_beta * (((_disp[_qp] - _disp_old[_qp])/(_dt*_dt)) - _vel_old[_qp]/_dt - accel_old * ( 0.5 - _beta));
}
