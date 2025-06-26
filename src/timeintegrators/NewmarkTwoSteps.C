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
/****************************************************************/

#include "NewmarkTwoSteps.h"
#include "NonlinearSystem.h"
#include "FEProblem.h"

registerMooseObject("MECHApp", NewmarkTwoSteps);

template<>
InputParameters validParams<NewmarkTwoSteps>()
{
  InputParameters params = validParams<TimeIntegrator>();
    params.addRequiredParam<Real>("rho_f", "density_fluid");
    params.addRequiredParam<Real>("rho_s", "density_solid");

  return params;
}

NewmarkTwoSteps::NewmarkTwoSteps(const InputParameters & parameters) :
    TimeIntegrator(parameters),
    _residual_old(_nl.addVector("residual_old", false, GHOSTED)),
    _residual_older(_nl.addVector("residual_older", false, GHOSTED)),
    _rho_s(getParam<Real>("rho_s")),
    _rho_f(getParam<Real>("rho_f")),    
    _solution_older(_sys.solutionState(2))
{
}


void
NewmarkTwoSteps::computeTimeDerivatives()
{
    // This is actually the second derivative, let us hope it works...

   NumericVector<Number> & _u_dot = *_sys.solutionUDot();
    
  _u_dot  = *_solution;

  //computeTimeDerivativeHelper(_u_dot, _solution_old);

  _u_dot -= _solution_old;
  _u_dot -= _solution_old;

  _u_dot += _solution_older;
  _u_dot *= (_rho_s - _rho_f) * 4.0/ _dt/_dt;
  _u_dot.close();

  _du_dot_du = (_rho_s - _rho_f) * 4.0/ _dt/_dt;
}

void
NewmarkTwoSteps::init()
{
  //if (_t_step == 1)
  {
    // make sure that time derivative contribution is zero in the first pre-solve step
    NumericVector<Number> & _u_dot = *_sys.solutionUDot();
    _u_dot.zero();
    _u_dot.close();

    _du_dot_du = (_rho_s -_rho_f) * 4.0/_dt/_dt;


    // for the first time step, compute residual for the old time step
    _nl.computeResidualTag(_nl.RHS(), _nl.nonTimeVectorTag());
    _residual_old = _nl.RHS();
    //_residual_old.close();
  }
}

void
NewmarkTwoSteps::postResidual(NumericVector<Number> & residual)
{
  if (_t_step == 1)
  {

      residual += _Re_non_time;
      residual += _residual_old;
      residual += _Re_time;
  }
    else
    {
        
        residual += _Re_non_time;
        residual += _residual_old;
        residual += _residual_old;
        residual += _residual_older;
        residual += _Re_time;
    }
}

void
NewmarkTwoSteps::postStep()
{
  // shift the residual in time
  _residual_older =_residual_old;
  //_residual_older.close();
  _residual_old = _Re_non_time;
  //_residual_old.close();

}

void
NewmarkTwoSteps::computeADTimeDerivatives(DualReal & ad_u_dot, const dof_id_type & dof) const
{
}

