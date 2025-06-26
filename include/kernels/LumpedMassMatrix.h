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

#ifndef LumpedMassMatrix_H
#define LumpedMassMatrix_H

#include "Kernel.h"

// Forward Declaration
class LumpedMassMatrix;

template<>
InputParameters validParams<LumpedMassMatrix>();

class LumpedMassMatrix : public Kernel
{
public:
  LumpedMassMatrix(const InputParameters & parameters);

protected:

  const VariableValue & _u_nodal;
    
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual void computeJacobian() override;
  virtual void computeResidual() override;
    

};
#endif 
