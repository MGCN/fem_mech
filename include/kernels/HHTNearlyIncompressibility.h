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
/*        MECH - ICS Imm Bound  simulation framework            */
/*                Prepared by Maria Nestola                     */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/*         Kernel for the Nearly Incompressibility              */
/****************************************************************/

#ifndef HHTNearlyIncompressibility_H
#define HHTNearlyIncompressibility_H

#include "Kernel.h"
#include "ElasticMaterial.h"

class HHTNearlyIncompressibility;

template <>
InputParameters validParams<HHTNearlyIncompressibility>();

class HHTNearlyIncompressibility : public Kernel
{
public:
  HHTNearlyIncompressibility(InputParameters const & parameters);

  unsigned int component;

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  unsigned const _numComp;
  
  unsigned int _ndisp;
  std::vector< unsigned int> _disp_var;
    
  RealVectorValue zeros;

  RealTensorValue H;
    
  MaterialProperty<Real> const &_J;
  MaterialProperty<Real> const &_J_old;
  MaterialProperty<RealTensorValue> const &_F;
  MaterialProperty<RealTensorValue> const &_invFtr;
  MaterialProperty<RealTensorValue> const &_invFtr_old;
  Real _kappa;
  Real _alpha_f;
  RealTensorValue _Cinv;
     
  AuxVariableName _var_name;
        
  Real ** _J_lin;
  RealTensorValue **_Pvol_Lin;
    
  unsigned int _jvar_index;
    
};

#endif // STRESSDIVERGENCELIBMESHTENSOR_H
