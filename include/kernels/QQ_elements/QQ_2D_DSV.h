#ifndef QQ_2D_DSV_H
#define QQ_2D_DSV_H

#include "Kernel.h"

class QQ_2D_DSV;

template<>
InputParameters validParams<QQ_2D_DSV>();


class QQ_2D_DSV : public Kernel
{
public:
    QQ_2D_DSV( InputParameters const & params);
protected:
    virtual Real computeQpResidual(){return 0.0;};
    //virtual Real computeQpJacobian(){return 0.0;};
    //virtual Real computeQpOffDiagJacobian(unsigned int jvar){std::cout<<jvar<<std::endl; return 0.0;};

    virtual void computeResidual();
    //virtual void computeJacobian(){};
    //virtual void computeOffDiagJacobian(unsigned int jvar);
    //virtual void computeOffDiagJacobianScalar(unsigned int jvar){std::cout<<jvar<<std::endl;};
    
    void computeResidual2D();
    //void computeJacobian2D();
    
    void computeGradient(RealVectorValue x0, RealVectorValue x1, RealVectorValue x2, RealVectorValue * Gradient);

    int _component;
    
    // indices of the displacement variables
    unsigned int _disp_x_var;
    unsigned int _disp_y_var;
//    unsigned int _disp_z_var;

    // Reference to local values of the displacement variables
    // In order to access to the right values, this kernel relies on the fact
    // that we employ Simposon quadrature rule.
    VariableValue    const & _disp_x;
    VariableValue    const & _disp_y;
//    VariableValue    const & _disp_z;

    RealTensorValue _identity;
    Real _mu,_lambda;

    // Here we store the map between simposn and P2 elements
    int * simpson_to_tri6;
    
    // Here we have the values of displacemnts in the right order
    Real * _u_x;
    Real * _u_y;
    
    // Local to global map
    // First index refers to the triangle 0:3,
    // second index refers to the node 0:2
    int **_local_to_global;
    
    // Here we store the gradients of the P1 functions
    // Using second order element, they are not available in the kernel
    // We employ Drosos's secret procedure
    RealVectorValue **_qq_grad_test;
    
    // Displacement gradient fine/coarse
    RealTensorValue *U;
    // Deformation gradient fine/coarse
    RealTensorValue *F;
    // Cauchy tensot fine/coarse
    RealTensorValue *C;
    // Strain fine/coarse
    RealTensorValue *E;
    // Strain quasi quadratic
    RealTensorValue *EQQ;
    // Stress quasi quadratic
    RealTensorValue *SQQ;
    
    Real * _v_x;
    Real * _v_y;
    
    RealTensorValue * _V;
    RealTensorValue * _E_lin;
    RealTensorValue * _E_lin_QQ;
};

#endif 
