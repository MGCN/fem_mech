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

#ifndef NEWMARKVELAUX_H
#define NEWMARKVELAUX_H

#include "AuxKernel.h"

class NewmarkVelAux;

template<>
InputParameters validParams<NewmarkVelAux>();

class NewmarkVelAux : public AuxKernel
{
public:
    
    /**
     *Calcualtes velocity using Newmark time integration scheme
     */
    NewmarkVelAux(const InputParameters & parameters);
    
    virtual ~NewmarkVelAux() {}
    
protected:
    virtual Real computeValue();
    
    
    const VariableValue & _accel_old;
    const VariableValue & _accel;
    
    Real _gamma;

     const VariableValue & u_old;
    
};

#endif //NEWMARKVELAUX_H
