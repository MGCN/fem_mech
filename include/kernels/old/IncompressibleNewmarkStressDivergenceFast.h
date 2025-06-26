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
/*                Prepared by Marco Favino,                     */
/*                  ICS, USI, 6900 Lugano                       */
/****************************************************************/



#ifndef IncompressibleNewmarkStressDivergenceFast_H
#define IncompressibleNewmarkStressDivergenceFast_H

#include "Kernel.h"
#include "ElasticMaterial.h"

class IncompressibleNewmarkStressDivergenceFast;

template <>
InputParameters validParams<IncompressibleNewmarkStressDivergenceFast>();

class IncompressibleNewmarkStressDivergenceFast : public Kernel
{
public:
  IncompressibleNewmarkStressDivergenceFast(
                                InputParameters const & parameters);



protected:

    virtual void computeResidual();
    virtual void computeResidual2D();
    virtual void computeResidual3D();
    
    virtual  void computeJacobian(){};
    virtual  void computeJacobian2D();
    virtual  void computeJacobian3D();
    
    virtual  void computeOffDiagJacobian(unsigned int );
    
    virtual Real computeQpResidual(){return 0.0;};

    virtual Real computeQpJacobian(){return 0.0;};

    virtual Real computeQpOffDiagJacobian(unsigned int){return 0.0;};
    
    
    unsigned const _dim;

    unsigned const _disp_x_var;
    unsigned const _disp_y_var;
    unsigned const _disp_z_var;
    unsigned const _pres_var;
    
    RealVectorValue zeros;

    RealTensorValue H;

    const VariableValue & _p;
    const VariableValue & _p_old;
    const VariableValue & _p_older;

    MaterialProperty<Real> const &_J;
    MaterialProperty<Real> const &_J_old;
    MaterialProperty<Real> const &_J_older;
    
    MaterialProperty<RealTensorValue> const &_invFtr;
    const MaterialProperty<RealTensorValue> &_invFtr_old;
    const MaterialProperty<RealTensorValue> &_invFtr_older;
    
    const MaterialProperty<RealTensorValue> &_P;
    const MaterialProperty<RealTensorValue> &_P_old;
    const MaterialProperty<RealTensorValue> &_P_older;
    
    RealTensorValue ***_stress_lin;

    
    const MaterialProperty<ElasticMaterial *> &_materialPointer;
    
    MooseVariable & _variable_pres;
    
    /// the current shape functions
    const VariableTestValue & _test_pres;

  
};

#endif // IncompressibleNewmarkStressDivergenceFast_H
