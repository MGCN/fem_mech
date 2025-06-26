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

#ifndef MASSBALANCEONLYRESIDUAL_H
#define MASSBALANCEONLYRESIDUAL_H

#include "Kernel.h"
#include "ElasticMaterial.h"

class MassBalanceOnlyResidual;

template <>
InputParameters validParams<MassBalanceOnlyResidual>();

class MassBalanceOnlyResidual : public Kernel
{
public:
  MassBalanceOnlyResidual(
                                InputParameters const & parameters);

  unsigned int _component;

protected:
  virtual Real computeQpResidual();

    virtual Real computeQpJacobian(){return 0.0;};

    virtual Real computeQpOffDiagJacobian(unsigned int jvar){return 0.0;};

    virtual void computeJacobian();
    
    virtual void computeOffDiagJacobian(unsigned int );
    
    unsigned const _numComp;
    

    MaterialProperty<Real> const &_J;
    
};

#endif // STRESSDIVERGENCELIBMESHTENSOR_H
