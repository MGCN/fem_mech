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

#ifndef AUXGROWTHFACTOR_H
#define AUXGROWTHFACTOR_H

#include "AuxKernel.h"

class AuxGrowthFactor;

template<>
InputParameters validParams<AuxGrowthFactor>();

class AuxGrowthFactor : public AuxKernel
{
public:
    
    /**
     *Calcualtes velocity using Newmark time integration scheme
     */
    AuxGrowthFactor(const InputParameters & parameters);
    
    virtual ~AuxGrowthFactor() {}
    
protected:
    virtual Real computeValue();
    
    
    // const VariableValue & _u_old;
   
    
      //const VariableValue & _u_old;
	Real _k_plus;
	Real _k_minus;

	const VariableValue & _threshold;

	Real _theta_plus;
	Real _theta_minus;
	Real _m_plus;
	Real _m_minus;
   const VariableValue & u_old;
    
};

#endif //NEWMARKVELAUX_H
