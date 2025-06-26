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
/*               Prepared by Maria Nestola,                     */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/****************************************************************/



#ifndef SplitNewmarkStressDivergenceFast_H
#define SplitNewmarkStressDivergenceFast_H

#include "Kernel.h"
#include "ElasticMaterial.h"

class SplitNewmarkStressDivergenceFast;

template <>
InputParameters validParams<SplitNewmarkStressDivergenceFast>();

class SplitNewmarkStressDivergenceFast : public Kernel
{
public:
    SplitNewmarkStressDivergenceFast(
                                        InputParameters const & parameters);
    
    
    
protected:
    
    virtual void computeResidual();
    
    virtual  void computeJacobian();
    
    virtual  void computeOffDiagJacobian(unsigned int );
    
    virtual Real computeQpResidual(){return 0.0;};
    
    virtual Real computeQpJacobian(){return 0.0;};
    
    virtual Real computeQpOffDiagJacobian(unsigned int){return 0.0;};
    
    unsigned const _dim;
    
    unsigned const _disp_x_var;
    unsigned const _disp_y_var;
    unsigned const _disp_z_var;
    unsigned * _disp_var;
    unsigned const _pres_var;
    
    RealVectorValue zeros;
    
    RealTensorValue H;
    
    bool const _has_pressure;
    const VariableValue & _p;
    
    
    MaterialProperty<Real> const &_J;
    MaterialProperty<RealTensorValue> const &_invFtr;
    
    const MaterialProperty<RealTensorValue> &_P;
    const MaterialProperty<RealTensorValue> &_P_old;
    
    RealTensorValue ***_stress_lin;
    
    
    const MaterialProperty<ElasticMaterial *> &_materialPointer;

    // const VariableValue & _reaction_x;
    
    // const VariableValue & _reaction_y;
    
    // const VariableValue & _reaction_z;
    
    // bool _contact;

    MooseVariable * _variable_pres;
    
    /// the current shape functions
    const VariableTestValue * _test_pres;
    
    DenseVector<Number> ** _residual;
    
    DenseVector<Number> * _local_re;
    
    DenseMatrix<Number> *** _jacobian;
    
    DenseMatrix <Number> ** _local_A;
    
    DenseMatrix<Number> ** _BT;
    
    DenseMatrix <Number> * _local_BT;
    
    DenseMatrix<Number> ** _B;
    
    DenseMatrix <Number> * _local_B;
    
    // RealVectorValue _reaction;
    
};

#endif // NewmarkStressDivergenceFast_H
