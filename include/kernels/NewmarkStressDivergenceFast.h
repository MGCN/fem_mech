///****************************************************************/
///*               DO NOT MODIFY THIS HEADER                      */
///* MOOSE - Multiphysics Object Oriented Simulation Environment  */
///*                                                              */
///*           (c) 2010 Battelle Energy Alliance, LLC             */
///*                   ALL RIGHTS RESERVED                        */
///*                                                              */
///*          Prepared by Battelle Energy Alliance, LLC           */
///*            Under Contract No. DE-AC07-05ID14517              */
///*            With the U. S. Department of Energy               */
///*                                                              */
///*            See COPYRIGHT for full restrictions               */
///****************************************************************/
///****************************************************************/
///*               DO NOT MODIFY THIS HEADER                      */
///*       MECH - ICS Mechanical simulation framework             */
///*                Prepared by Marco Favino,                     */
///*                  ICS, USI, 6900 Lugano                       */
///*                                                              */
///* Fast Kernel for the Incompressibility Stress Divergence.     */
///* To be used with the EmptyKernel.                             */
///****************************************************************/
//
//
//
//#ifndef NewmarkStressDivergenceFast_H
//#define NewmarkStressDivergenceFast_H
//
//#include "Kernel.h"
//#include "ElasticMaterial.h"
//
//class NewmarkStressDivergenceFast;
//
//template <>
//InputParameters validParams<NewmarkStressDivergenceFast>();
//
//class NewmarkStressDivergenceFast : public Kernel
//{
//public:
//  NewmarkStressDivergenceFast(
//                                InputParameters const & parameters);
//
//
//
//protected:
//
//    virtual void computeResidual();
//    virtual void computeResidual2D();
//    virtual void computeResidual3D();
//
//    virtual  void computeJacobian(){};
//    virtual  void computeJacobian2D();
//    virtual  void computeJacobian3D();
//
//    virtual  void computeOffDiagJacobian(unsigned int );
//
//    virtual Real computeQpResidual(){return 0.0;};
//
//    virtual Real computeQpJacobian(){return 0.0;};
//
//    virtual Real computeQpOffDiagJacobian(unsigned int){return 0.0;};
//
//
//    unsigned const _dim;
//
//    unsigned const _disp_x_var;
//    unsigned const _disp_y_var;
//    unsigned const _disp_z_var;
//    unsigned const _pres_var;
//
//    RealVectorValue zeros;
//
//    RealTensorValue H;
//
//    const VariableValue & _p;
//    const VariableValue & _p_old;
//    const VariableValue & _p_older;
//
//    MaterialProperty<Real> const &_J;
//    MaterialProperty<Real> const &_J_old;
//    MaterialProperty<Real> const &_J_older;
//
//    MaterialProperty<RealTensorValue> const &_invFtr;
//    const MaterialProperty<RealTensorValue> &_invFtr_old;
//    const MaterialProperty<RealTensorValue> &_invFtr_older;
//
//    const MaterialProperty<RealTensorValue> &_P;
//    const MaterialProperty<RealTensorValue> &_P_old;
//    const MaterialProperty<RealTensorValue> &_P_older;
//
//    RealTensorValue ***_stress_lin;
//
//
//    const MaterialProperty<ElasticMaterial *> &_materialPointer;
//
//    MooseVariable * _variable_pres;
//
//
//    /// the current shape functions
//    const VariableTestValue *_test_pres;
//
//
//
//};
//
//#endif // NewmarkStressDivergenceFast_H
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



#ifndef NewmarkSTRESSDIVERGENCEFAST_H
#define NewmarkSTRESSDIVERGENCEFAST_H

#include "Kernel.h"
#include "ElasticMaterial.h"

class NewmarkStressDivergenceFast;

template <>
InputParameters validParams<NewmarkStressDivergenceFast>();

class NewmarkStressDivergenceFast : public Kernel
{
public:
    NewmarkStressDivergenceFast(InputParameters const & parameters);
    ~NewmarkStressDivergenceFast();
    
    
    
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
    unsigned  _pres_var;
    
    RealVectorValue zeros;
    
    RealTensorValue H;
    
    bool const _has_pressure;
    bool _is_linear;
    MaterialProperty<Real> const & _p;
    MaterialProperty<Real> const & _J;
    MaterialProperty<RealTensorValue> const &_invFtr;
    
   const MaterialProperty<RealTensorValue> &_P;
   const MaterialProperty<RealTensorValue> &_P_old;
   const MaterialProperty<RealTensorValue> &_P_older;
    
    RealTensorValue ***_stress_lin;
    
    
    const MaterialProperty<ElasticMaterial *> &_materialPointer;
    MaterialProperty<RealTensorValue> const &_gradDisp;
    
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
    
    
    
};

#endif // STRESSDIVERGENCEFAST_H
