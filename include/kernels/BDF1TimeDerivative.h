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

#ifndef BDF1TIMEDERIVATIVE_H
#define BDF1TIMEDERIVATIVE_H

#include "TimeDerivative.h"

// Forward Declarations
class BDF1TimeDerivative;

template<>
InputParameters validParams<BDF1TimeDerivative>();

class BDF1TimeDerivative : public TimeDerivative
{
public:

  BDF1TimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;

	const VariableValue & _u_old;

	const VariableValue & _u_older;

	Real & _dt;

	int & _t_step;
    
    Real _rho_s;
    
    Real _rho_f;
};

#endif //BDF1TIMEDERIVATIVE_H
