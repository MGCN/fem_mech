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

#ifndef EULERVELAUX_H
#define EULERVELAUX_H

#include "AuxKernel.h"

class EulerVelAux;

template<>
InputParameters validParams<EulerVelAux>();

class EulerVelAux : public AuxKernel
{
public:
    
    /**
     *Calcualtes velocity using Newmark time integration scheme
     */
    EulerVelAux(const InputParameters & parameters);
    
    virtual ~EulerVelAux() {}
    
protected:
    virtual Real computeValue();
    
    
    const VariableValue & _u_old;
    const VariableValue & _u;
    
    Real _gamma;
    
};

#endif //NEWMARKVELAUX_H
