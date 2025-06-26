//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef CONSERVATIVEADVECTIONSUPG_H
#define CONSERVATIVEADVECTIONSUPG_H

#include "Kernel.h"

// Forward Declaration
class ConservativeAdvectionSUPG;

/**
 * Advection of the variable by the velocity provided by the user.
 * Options for numerical stabilization are: none; full upwinding
 */
template <>
InputParameters validParams<ConservativeAdvectionSUPG>();

class ConservativeAdvectionSUPG : public Kernel
{
public:
    ConservativeAdvectionSUPG(const InputParameters & parameters);
    
protected:
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;
    /// advection velocity
    //RealVectorValue _velocity;
    const VariableGradient &_vel;
    Real _k;
    
    Real negSpeedQp() const;
    

};

#endif // CONSERVATIVEADVECTION_H
