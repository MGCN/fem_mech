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

#ifndef NEWMARKTWOSTEPS_H
#define NEWMARKTWOSTEPS_H

#include "TimeIntegrator.h"

class NewmarkTwoSteps;

template<>
InputParameters validParams<NewmarkTwoSteps>();

/**
 * Crank-Nicolson time integrator.
 *
 * The scheme is defined as:
 *   \f$ \frac{du}{dt} = 1/2 * (F(U^{n+1}) + F(U^{n})) \f$,
 * but the form we are using it in is:
 *   \f$ 2 * \frac{du}{dt} = (F(U^{n+1}) + F(U^{n})) \f$.
 */
class NewmarkTwoSteps : public TimeIntegrator
{
public:
  NewmarkTwoSteps(const InputParameters & parameters);
    virtual ~NewmarkTwoSteps(){};

    virtual void init() override;
    virtual int order() override { return 2; }
    virtual void computeTimeDerivatives() override;
    virtual void postResidual(NumericVector<Number> & residual) override;
    virtual void postStep() override;
    virtual void computeADTimeDerivatives(DualReal & ad_u_dot, const dof_id_type & dof) const override ;

protected:
  NumericVector<Number> & _residual_old;
  NumericVector<Number> & _residual_older;
    
    Real _rho_s;
    Real _rho_f;
  const NumericVector<Number> & _solution_older;
};

#endif /* CRANKNICOLSON_H_ */
