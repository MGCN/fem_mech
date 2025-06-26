//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html



#ifndef INERTIALFORCEMECH_H
#define INERTIALFORCEMECH_H



#include "Kernel.h"



// Forward Declarations
class InertialForceMech;



template <>
InputParameters validParams<InertialForceMech>();



class InertialForceMech : public Kernel
{
public:
InertialForceMech(const InputParameters & parameters);



protected:
virtual Real computeQpResidual();



virtual Real computeQpJacobian();



private:
const Real _density;
const VariableValue & _u_old;
const VariableValue & _u_older;
const Real _eta;
};



#endif // InertialForceMech_H
