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
/*                                                              */
/****************************************************************/



#ifndef CONVECTIONNS_H
#define CONVECTIONNS_H

#include "Kernel.h"
#include "ElasticMaterial.h"

class ConvectionNS;

template <>
InputParameters validParams<ConvectionNS>();

class ConvectionNS : public Kernel
{
public:
 ConvectionNS(InputParameters const & parameters);



protected:

    virtual void computeResidual();
    
    virtual  void computeJacobian();
    
    virtual  void computeOffDiagJacobian(unsigned int );
    
    virtual Real computeQpResidual(){return 0.0;};

    virtual Real computeQpJacobian(){return 0.0;};

    virtual Real computeQpOffDiagJacobian(unsigned int){return 0.0;};
    
    unsigned const _dim;
    
    unsigned const _vel_x_var;
    unsigned const _vel_y_var;
    unsigned const _vel_z_var;

    MaterialProperty<RealTensorValue> const &_grad_vel;
    VariableValue const & _vel_y;
    VariableValue const & _vel_z;
    
    unsigned * _vel_var;
    
    RealTensorValue _H;
    RealVectorValue _h;
    
    
    RealVectorValue *_vel;
    
    RealVectorValue *_non_linear_term;
    
    RealVectorValue ***_stress_lin;
    
    DenseVector<Number> ** _residual;
    
    DenseVector<Number> * _local_re;

    DenseMatrix<Number> *** _jacobian;
    
    DenseMatrix <Number> ** _local_A;

};

#endif // STRESSDIVERGENCEFAST_H
