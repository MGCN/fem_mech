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
/* Fast Kernel for the Incompressibility Stress Divergence.     */
/* To be used with the EmptyKernel.                             */
/****************************************************************/



#include "CompressibleNewmarkStressDivergenceFast.h"

#include "Material.h"
#include "Assembly.h"
#include "MooseMesh.h"

#include "libmesh/quadrature.h"

registerMooseObject("MECHApp", CompressibleNewmarkStressDivergenceFast);

template <>
InputParameters validParams<CompressibleNewmarkStressDivergenceFast>() {

  // inherit the parameters of the Kernels:
  InputParameters params = validParams<Kernel>();

  params.addRequiredCoupledVar("disp_x","disp_x");
  params.addRequiredCoupledVar("disp_y","disp_y");
  params.addCoupledVar("disp_z","disp_z");
  return params;
}


CompressibleNewmarkStressDivergenceFast::CompressibleNewmarkStressDivergenceFast(
                                         InputParameters const & parameters) :

    // inherit the parameters of the Kernels:
    Kernel(parameters),
    // We are doing elasticity from R^n to R^n,
    // hence the number of components is the same as the mesh dimension:
    // For example for _mesh.dimension=3 we have
    //_component = 0 for disp_x, _component = 1 for disp_y, _component = 2 for
    // disp_z
    _dim(_mesh.dimension()),
    
    // the pressure variable in coupled and non coupled copy:
    _disp_x_var(coupled("disp_x")),
    _disp_y_var(coupled("disp_y")),
    _disp_z_var(_mesh.dimension() == 3 ? coupled("disp_z") : 100000),
    _P(getMaterialProperty<RealTensorValue>("firstPiolaKirchhofTensor")),
    _P_old(getMaterialPropertyOld<RealTensorValue>("firstPiolaKirchhofTensor")),
    _P_older(getMaterialPropertyOlder<RealTensorValue>("firstPiolaKirchhofTensor")),
    _materialPointer(getMaterialProperty<ElasticMaterial *>("materialPointer"))


{

  // we suppose that the components of displacements are the first three.
  if (_dim == 1)
  {
     mooseError( "Stress divergence fast does not work with dim=1");
  }

    for (int i = 0; i < 3; ++i)
        zeros(i) = 0.0;
    
    
    _stress_lin = new RealTensorValue **[_dim] /*[_phi.size()]*/;
    
    for (int d=0; d<_dim; ++d)
        _stress_lin[d]=new RealTensorValue *[64];
    
    
    for (int d=0; d<_dim; ++d)
        for (int j = 0; j < 64 /*_phi.size()*/; ++j)
            _stress_lin[d][j] = new RealTensorValue [64] /*[_qrule->n_points()]*/;

    
}

void
CompressibleNewmarkStressDivergenceFast::computeOffDiagJacobian(unsigned int jvar)
{
    //DenseMatrix<Number> & ke00 = _assembly.jacobianBlock( jvar);
    if (jvar==_disp_x_var)
    {
        if (_dim==2)
            computeJacobian2D();
        else
        {
            if (_dim==3)
                computeJacobian3D();
            else
                mooseError( "You should not be here");
        }
    }
}


void
CompressibleNewmarkStressDivergenceFast::computeJacobian2D()
{
    DenseMatrix<Number> & A00 = _assembly.jacobianBlock(_disp_x_var, _disp_x_var);
    DenseMatrix<Number> & A01 = _assembly.jacobianBlock(_disp_x_var, _disp_y_var);

    
    DenseMatrix<Number> & A10 = _assembly.jacobianBlock(_disp_y_var, _disp_x_var);
    DenseMatrix<Number> & A11 = _assembly.jacobianBlock(_disp_y_var, _disp_y_var);

    
    
    DenseMatrix<Number> _local_A[2][2];

    
    // set dimensions of local Jacobians and set them to zero
    for (int i=0; i<2; ++i)
        for (int j=0; j<2; ++j)
        {
            _local_A[i][j].resize(A00.m(),A00.n());
            _local_A[i][j].zero();
        }


    
    
    for (int dim=0; dim<2; ++dim)
    {
        H = RealTensorValue(zeros, zeros, zeros);
        for (int j=0; j<_phi.size(); ++j)
        {
            for (int qp=0; qp<_qrule->n_points(); ++qp)
            {
                for (int k=0 ; k<2; ++k)
                    H(dim,k)=_grad_phi[j][qp](k);
                (*(_materialPointer[qp])).evaluate_stress_lin(qp, H, _stress_lin[dim][j][qp]);
                _stress_lin[dim][j][qp] = 1./4 *_stress_lin[dim][j][qp];
            }
        }
    }
    
    for (int idim=0; idim<2; ++idim)
        for (int jdim=idim; jdim<2; ++jdim)
            for (_i = 0; _i < _test.size(); _i++)
                for (_j = 0; _j < _phi.size(); _j++)
                    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
                        _local_A[idim][jdim](_i,_j)+=_JxW[_qp] * _coord[_qp]*
                        (_stress_lin[jdim][_j][_qp](idim,0)*_grad_test[_i][_qp](0)+
                         _stress_lin[jdim][_j][_qp](idim,1)*_grad_test[_i][_qp](1)  // forse jdim?
                         );
    

    
    
    _local_A[0][1].get_transpose(_local_A[1][0]);

    A00 += _local_A[0][0];
    A01 += _local_A[0][1];
    A10 += _local_A[1][0];
    A11 += _local_A[1][1];

    
    
    if(_has_diag_save_in)
    {
        mooseError("Error: diag in already saved in");
    }
    
    
}

void
CompressibleNewmarkStressDivergenceFast::computeJacobian3D()
{
//    DenseMatrix<Number> & ke00 = _assembly.jacobianBlock(0, 0);
//    DenseMatrix<Number> & ke01 = _assembly.jacobianBlock(0, 1);
//    DenseMatrix<Number> & ke02 = _assembly.jacobianBlock(0, 2);
//    
//    DenseMatrix<Number> & ke10 = _assembly.jacobianBlock(1, 0);
//    DenseMatrix<Number> & ke11 = _assembly.jacobianBlock(1, 1);
//    DenseMatrix<Number> & ke12 = _assembly.jacobianBlock(1, 2);
//
//    DenseMatrix<Number> & ke20 = _assembly.jacobianBlock(2, 0);
//    DenseMatrix<Number> & ke21 = _assembly.jacobianBlock(2, 1);
//    DenseMatrix<Number> & ke22 = _assembly.jacobianBlock(2, 2);
//
//    
//    DenseMatrix<Number> _local_ke[3][3];
//    
//    for (int i=0; i<3; ++i)
//        for (int j=0; j<3; ++j)
//        {
//            _local_ke[i][j].resize(ke00.m(),ke00.n());
//            _local_ke[i][j].zero();
//        }
//    
//    
//    for (int dim=0; dim<3; ++dim)
//    {
//        H = RealTensorValue(zeros, zeros, zeros);
//        for (int j=0; j<_phi.size(); ++j)
//        {
//            for (int qp=0; qp<_qrule->n_points(); ++qp)
//            {
//                for (int k=0 ; k<3; ++k)
//                    H(dim,k)=_grad_phi[j][qp](k);
//                (*(_materialPointer[qp])).evaluate_stress_lin(qp, H, _stress_lin[dim][j][qp]);
//                _stress_lin[dim][j][qp]+=_p[qp]*_J[qp]* ( _invFtr[qp].contract(H)*_invFtr[qp]-_invFtr[qp]*H.transpose()*_invFtr[qp]  );
//            }
//        }
//    }
//    
//    for (int idim=0; idim<3; ++idim)
//        for (int jdim=idim; jdim<3; ++jdim)
//            for (_i = 0; _i < _test.size(); _i++)
//                for (_j = 0; _j < _phi.size(); _j++)
//                    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
//                        _local_ke[idim][jdim](_i,_j)+=_JxW[_qp] * _coord[_qp]*
//                        (_stress_lin[jdim][_j][_qp](idim,0)*_grad_test[_i][_qp](0)+
//                         _stress_lin[jdim][_j][_qp](idim,1)*_grad_test[_i][_qp](1)+
//                         _stress_lin[jdim][_j][_qp](idim,2)*_grad_test[_i][_qp](2)  // forse jdim?
//                         );
//    
//    
//    for (int i=0; i<_local_ke[1][0].m(); ++i)
//        for (int j=0; j<_local_ke[1][0].m(); ++j)
//            _local_ke[1][0](i,j)=_local_ke[0][1](j,i);
//
//    for (int i=0; i<_local_ke[1][0].m(); ++i)
//        for (int j=0; j<_local_ke[1][0].m(); ++j)
//            _local_ke[2][0](i,j)=_local_ke[0][2](j,i);
//
//    for (int i=0; i<_local_ke[1][0].m(); ++i)
//        for (int j=0; j<_local_ke[1][0].m(); ++j)
//            _local_ke[2][1](i,j)=_local_ke[1][2](j,i);
//
//    
//    ke00 += _local_ke[0][0];
//    ke01 += _local_ke[0][1];
//    ke02 += _local_ke[0][2];
//    
//    ke10 += _local_ke[1][0];
//    ke11 += _local_ke[1][1];
//    ke12 += _local_ke[1][2];
//    
//    ke20 += _local_ke[2][0];
//    ke21 += _local_ke[2][1];
//    ke22 += _local_ke[2][2];
//
//
//    
//    if(_has_diag_save_in)
//    {
//        mooseError("Error: diag in already saved in");
//    }

    
}


void CompressibleNewmarkStressDivergenceFast::computeResidual()
{
    if (_dim==2)
        computeResidual2D();
    else
    {
        if (_dim==3)
            computeResidual3D();
        else
            mooseError( "You should not be here");
        
    }
}

void
CompressibleNewmarkStressDivergenceFast::computeResidual2D()
{
    // References to the residual components
    DenseVector<Number> & re_x = _assembly.residualBlock(_disp_x_var);
    DenseVector<Number> & re_y = _assembly.residualBlock(_disp_y_var);
    
    // Local vectors where the local contribution is assembled
    DenseVector<Number> _local_re[2];
    for (int i=0; i<2; ++i)
    {
        _local_re[i].resize(re_x.size());
        _local_re[i].zero();
    }

    
    // Assembly of residual for mechanics
    for (int d=0; d<2; ++d)
        for (_i = 0; _i < _test.size(); ++_i)
            for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
            {
                if (_t_step==1){
                    
                    RealTensorValue P= 1./4 * ( _P[_qp] +_P_old[_qp]);
                    
                    _local_re[d](_i) += _JxW[_qp] * _coord[_qp] *
                    (
                     P(d , 0) * _grad_test[_i][_qp] (0)+
                     P(d , 1) * _grad_test[_i][_qp] (1)
                     );
                }
                if (_t_step>1){
                    
                    RealTensorValue P= 1./4 * ( _P[_qp]) +
                                       1./2 * (_P_old[_qp])+
                                       1./4 * (_P_older[_qp]) ;
                    
                    _local_re[d](_i) += _JxW[_qp] * _coord[_qp] *
                                                                (
                                                                P(d , 0) * _grad_test[_i][_qp] (0)+
                                                                P(d , 1) * _grad_test[_i][_qp] (1)
                                                                );
                }
            }
    
    

    
    // sum the local contributions to the global residual
    re_x += _local_re[0];
    re_y += _local_re[1];

}

void
CompressibleNewmarkStressDivergenceFast::computeResidual3D()
{
//    DenseVector<Number> & re_x = _assembly.residualBlock(_disp_x_var);
//    DenseVector<Number> & re_y = _assembly.residualBlock(_disp_y_var);
//    DenseVector<Number> & re_z = _assembly.residualBlock(_disp_z_var);
//    DenseVector<Number> & re_p = _assembly.residualBlock(_pres_var);
//    
//    DenseVector<Number> _local_re[3];
//    for (int i=0; i<3; ++i)
//    {
//        _local_re[i].resize(re_x.size());
//        _local_re[i].zero();
//    }
//    DenseVector<Number> _local_re_p;
//    _local_re_p.resize(re_p.size());
//    
//    for (int d=0; d<3; ++d)
//        for (_i = 0; _i < _test.size(); _i++)
//            for (_qp = 0; _qp < _qrule->n_points(); _qp++)
//            {
//                RealTensorValue P=_P[_qp]+_p[_qp]*_J[_qp]*_invFtr[_qp];
//                
//                _local_re[d](_i) += _JxW[_qp] * _coord[_qp] *
//                (
//                 P(d , 0) * _grad_test[_i][_qp] (0)+
//                 P(d , 1) * _grad_test[_i][_qp] (1)+
//                 P(d , 2) * _grad_test[_i][_qp] (2)
//                 );
//            }
//    
//    re_x += _local_re[0];
//    re_y += _local_re[1];
//    re_z += _local_re[2];
//    re_p += _local_re_p;
//    
}
