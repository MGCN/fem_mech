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
/*        MECH - ICS Imm Bound  simulation framework            */
/*                Prepared by Maria Nestola                     */
/*                                                              */	
/*      Kernel to consider the intertial term.                  */
/*      It implements the term                                  */
/*     _int ((_rho_s - _rho_f) * u_tt * test * dV )             */
/*      by using the Newmark method                             */
/****************************************************************/

#include "SplitNewmarkVelocity.h"
#include "SubProblem.h"
#include "MooseMesh.h"

registerMooseObject("MECHApp", SplitNewmarkVelocity);

template<>
InputParameters validParams<SplitNewmarkVelocity>()
{
  InputParameters params = validParams<Kernel>();
  // params.addRequiredParam<Real>("rho_f", "density_fluid");
  // params.addRequiredParam<Real>("rho_s", "density_solid");
  params.addRequiredCoupledVar("rho","density");
  params.addRequiredParam<unsigned>("component", "component");
  // params.addRequiredCoupledVar("contact_force_x","contact_force_x");
  // params.addCoupledVar("contact_force_y","contact_force_y");
  // params.addCoupledVar("contact_force_z","contact_force_z");
  // params.addParam<bool>("contact",false,
  //                         "Set to false when you want to pure L2-projection.");

    

  return params;
}

SplitNewmarkVelocity::SplitNewmarkVelocity(const InputParameters & parameters) :
    Kernel(parameters),
    _dim(_mesh.dimension()),
    // _rho_f(getParam<Real>("rho_f")),
    // _rho_s(getParam<Real>("rho_s")),
    _rho(coupledValue("rho")),
    _u_old(valueOld()),
    _component(getParam<unsigned>("component")),
    _P(getMaterialProperty<RealTensorValue>("firstPiolaKirchhofTensor")),
    _P_old(getMaterialPropertyOld<RealTensorValue>("firstPiolaKirchhofTensor"))
    // _reaction_x(coupledValueOld("contact_force_x")),
    // _reaction_y(_mesh.dimension() == 2 ? coupledValueOld("contact_force_y") : _zero),
    // _reaction_z(_mesh.dimension() == 3 ? coupledValueOld("contact_force_z") : _zero),
    // _contact(getParam<bool>("contact"))

{}

Real
SplitNewmarkVelocity::computeQpResidual()
{
    //Real diff_rho = _rho_s - _rho_f;
    
    P = 1./2 * ( _P[_qp]  + _P_old[_qp] );


   
  // RealVectorValue _reaction(0.0,0.0,0.0);
    
  // if(_dim==2)
  // _reaction=(_reaction_x[_qp],_reaction_y[_qp],0.0);
  // if(_dim==3)
  // _reaction=(_reaction_x[_qp],_reaction_y[_qp],_reaction_y[_qp]);
 
  //  if (_contact)
                           
   //  return _rho[_qp] * (_u[_qp] - _u_old[_qp])/_dt * _test[_i][_qp] +
   //                      P(_component , 0) * _grad_test[_i][_qp] (0) +
   //                      P(_component , 1) * _grad_test[_i][_qp] (1) +
   //                      P(_component , 2) * _grad_test[_i][_qp] (1) 
   //                      - _reaction(_component);
   // else
   
   return _rho[_qp] * (_u[_qp] - _u_old[_qp])/_dt * _test[_i][_qp] +
                      P(_component , 0) * _grad_test[_i][_qp] (0) +
                      P(_component , 1) * _grad_test[_i][_qp] (1) +
                      P(_component , 2) * _grad_test[_i][_qp] (1) ;

}

Real
SplitNewmarkVelocity::computeQpJacobian()
{
    //Real diff_rho =_rho_s - _rho_f;
    return 0.0;
    //diff_rho * _phi[_j][_qp]/_dt * _test[_i][_qp];

}
