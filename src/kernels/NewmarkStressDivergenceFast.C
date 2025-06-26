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
///*                Prepared by Maria Nestola,                    */
///*                  ICS, USI, 6900 Lugano                       */
///*                                                              */
///* Fast Kernel for the Incompressibility Stress Divergence.     */
///* To be used with the EmptyKernel.                             */
///****************************************************************/
//
//
//
//#include "NewmarkStressDivergenceFast.h"
//
//#include "Material.h"
//#include "Assembly.h"
//#include "MooseMesh.h"
//#include "MooseVariable.h"
//#include "libmesh/quadrature.h"
//
//
//registerMooseObject("MECHApp", NewmarkStressDivergenceFast);
//
//template <>
//InputParameters validParams<NewmarkStressDivergenceFast>() {
//
//  // inherit the parameters of the Kernels:
//  InputParameters params = validParams<Kernel>();
//
//  params.addRequiredCoupledVar("disp_x","disp_x");
//  params.addRequiredCoupledVar("disp_y","disp_y");
//  params.addCoupledVar("disp_z","disp_z");
//  // solve for the pressure in the massBalance Kernel implies pressure as coupled variable:
//  params.addCoupledVar("pressure","pressure");
//
//  params.addParam<NonlinearVariableName>("variable_pres", "Pressure variable");
//
//  return params;
//}
//
//
//NewmarkStressDivergenceFast::NewmarkStressDivergenceFast(
//                                         InputParameters const & parameters) :
//
//    // inherit the parameters of the Kernels:
//    Kernel(parameters),
//    // We are doing elasticity from R^n to R^n,
//    // hence the number of components is the same as the mesh dimension:
//    // For example for _mesh.dimension=3 we have
//    //_component = 0 for disp_x, _component = 1 for disp_y, _component = 2 for
//    // disp_z
//    _dim(_mesh.dimension()),
//
//    // the pressure variable in coupled and non coupled copy:
//    _disp_x_var(coupled("disp_x")),
//    _disp_y_var(coupled("disp_y")),
//    _disp_z_var(_mesh.dimension() == 3 ? coupled("disp_real_z") : 100000),
//    _p(getMaterialProperty<Real>("pressureFromMaterial")),
//    _p(parameters.isParamValid("pressure") ? coupledValue("pressure"):_zero),
//    _p_old(parameters.isParamValid("pressure") ? coupledValueOld("pressure"):_zero),
//    _p_older(parameters.isParamValid("pressure") ? coupledValueOlder("pressure"):_zero),
//    // inherit material properties:
//    _J(getMaterialProperty<Real>("deformationDeterminant")),
//    _J_old(getMaterialPropertyOld<Real>("deformationDeterminant")),
//    _J_older(getMaterialPropertyOlder<Real>("deformationDeterminant")),
//    _invFtr(getMaterialProperty<RealTensorValue>("invFtr")),
//    _invFtr_old(getMaterialPropertyOld<RealTensorValue>("invFtr")),
//    _invFtr_older(getMaterialPropertyOlder<RealTensorValue>("invFtr")),
//    _P(getMaterialProperty<RealTensorValue>("firstPiolaKirchhofTensor")),
//    _P_old(getMaterialPropertyOld<RealTensorValue>("firstPiolaKirchhofTensor")),
//    _P_older(getMaterialPropertyOlder<RealTensorValue>("firstPiolaKirchhofTensor")),
//    _materialPointer(getMaterialProperty<ElasticMaterial *>("materialPointer"))
//
//{
//    if(_has_pressure)
//    {
//        _variable_pres = &_subproblem.getStandardVariable(_tid, getParam<NonlinearVariableName>("pressure"));
//        _test_pres = &(_variable_pres[0].phi());
//        _pres_var=_variable_pres[0].number() ;
//    }
//
//
//
//  // we suppose that the components of displacements are the first three.
//  if (_dim == 1)
//  {
//     mooseError( "Stress divergence fast does not work with dim=1");
//  }
//
//    for (int i = 0; i < 3; ++i)
//        zeros(i) = 0.0;
//
//
//    _stress_lin = new RealTensorValue **[_dim] /*[_phi.size()]*/;
//
//    for (int d=0; d<_dim; ++d)
//        _stress_lin[d]=new RealTensorValue *[64];
//
//
//    for (int d=0; d<_dim; ++d)
//        for (int j = 0; j < 64 /*_phi.size()*/; ++j)
//            _stress_lin[d][j] = new RealTensorValue [64] /*[_qrule->n_points()]*/;
//
//
//}
//
//void
//NewmarkStressDivergenceFast::computeOffDiagJacobian(unsigned int jvar)
//{
//    //DenseMatrix<Number> & ke00 = _assembly.jacobianBlock( jvar);
//    if (jvar==_disp_x_var)
//    {
//        if (_dim==2)
//            computeJacobian2D();
//        else
//        {
//            if (_dim==3)
//                computeJacobian3D();
//            else
//                mooseError( "You should not be here");
//        }
//    }
//}
//
//
//void
//NewmarkStressDivergenceFast::computeJacobian2D()
//{
//    DenseMatrix<Number> & A00 = _assembly.jacobianBlock(_disp_x_var, _disp_x_var);
//    DenseMatrix<Number> & A01 = _assembly.jacobianBlock(_disp_x_var, _disp_y_var);
//
//
//    DenseMatrix<Number> & A10 = _assembly.jacobianBlock(_disp_y_var, _disp_x_var);
//    DenseMatrix<Number> & A11 = _assembly.jacobianBlock(_disp_y_var, _disp_y_var);
//
//
//
//    DenseMatrix<Number> _local_A[2][2];
//
//
//
//    DenseMatrix<Number> _local_BT[2];
//    DenseMatrix<Number> _local_B [2];
//
//    DenseMatrix<Number> & B12 = _assembly.jacobianBlock(_disp_y_var, _pres_var   );
//    DenseMatrix<Number> & B02 = _assembly.jacobianBlock(_disp_x_var, _pres_var  );
//    DenseMatrix<Number> & B20 = _assembly.jacobianBlock(_pres_var  , _disp_x_var);
//    DenseMatrix<Number> & B21 = _assembly.jacobianBlock(_pres_var  , _disp_y_var);
//
//
//
//    // set dimensions of local Jacobians and set them to zero
//    for (int i=0; i<2; ++i)
//        for (int j=0; j<2; ++j)
//        {
//            _local_A[i][j].resize(A00.m(),A00.n());
//            _local_A[i][j].zero();
//        }
//
//    if(_pres_var!= NULL){
//
//        for (int j=0; j<2; ++j)
//        {
//            _local_BT[j].resize(B02.m(),B02.n());
//            _local_BT[j].zero();
//        }
//
//
//        for (int j=0; j<2; ++j)
//        {
//            _local_B[j].resize(B20.m(),B20.n());
//            _local_B[j].zero();
//        }
//    }
//
//
//
//    for (int dim=0; dim<2; ++dim)
//    {
//        H = RealTensorValue(zeros, zeros, zeros);
//        for (int j=0; j<_phi.size(); ++j)
//        {
//            for (int qp=0; qp<_qrule->n_points(); ++qp)
//            {
//                for (int k=0 ; k<2; ++k)
//                    H(dim,k)=_grad_phi[j][qp](k);
//                (*(_materialPointer[qp])).evaluate_stress_lin(qp, H, _stress_lin[dim][j][qp]);
//                if(_pres_var!= NULL){_stress_lin[dim][j][qp]+=_p[qp]*_J[qp]* ( _invFtr[qp].contract(H)*_invFtr[qp]-_invFtr[qp]*H.transpose()*_invFtr[qp]);}
//                _stress_lin[dim][j][qp] = 1./4 *_stress_lin[dim][j][qp];
//            }
//        }
//    }
//
//    for (int idim=0; idim<2; ++idim)
//        for (int jdim=idim; jdim<2; ++jdim)
//            for (_i = 0; _i < _test.size(); _i++)
//                for (_j = 0; _j < _phi.size(); _j++)
//                    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
//                        _local_A[idim][jdim](_i,_j)+=_JxW[_qp] * _coord[_qp]*
//                        (_stress_lin[jdim][_j][_qp](idim,0)*_grad_test[_i][_qp](0)+
//                         _stress_lin[jdim][_j][_qp](idim,1)*_grad_test[_i][_qp](1)  // forse jdim?
//                         );
//
//
//    if(_pres_var!= NULL){
//        for (int idim=0; idim<2; ++idim)
//            for (_i=0; _i<_test.size(); ++_i)
//                for (_j=0; _j<_test_pres->size(); ++_j)
//                    for (_qp=0; _qp<_qrule->n_points(); ++_qp)
//                        _local_BT[idim](_i,_j)+= 1./4 * _JxW[_qp] * _coord[_qp]*((_test_pres[0][_j][_qp]))*_J[_qp]*
//                        (_invFtr[_qp](idim,0)*_grad_test[_i][_qp](0)+_invFtr[_qp](idim,1)*_grad_test[_i][_qp](1));
//    }
//
//
//    _local_A[0][1].get_transpose(_local_A[1][0]);
//
//    if(_pres_var!= NULL){
//        _local_BT[0].get_transpose(_local_B[0]);
//        _local_BT[1].get_transpose(_local_B[1]);
//    }
//
//    A00 += _local_A[0][0];
//    A01 += _local_A[0][1];
//    A10 += _local_A[1][0];
//    A11 += _local_A[1][1];
//
//    if(_pres_var!= NULL){
//
//        B02 += _local_BT[0];
//        B12 += _local_BT[1];
//
//        B20 += _local_B[0];
//        B21 += _local_B[1];
//    }
//
//
//    if(_has_diag_save_in)
//    {
//        mooseError("Error: diag in already saved in");
//    }
//
//
//}
//
//void
//NewmarkStressDivergenceFast::computeJacobian3D()
//{
//    DenseMatrix<Number> & A00 = _assembly.jacobianBlock(_disp_x_var, _disp_x_var);
//    DenseMatrix<Number> & A01 = _assembly.jacobianBlock(_disp_x_var, _disp_y_var);
//    DenseMatrix<Number> & A02 = _assembly.jacobianBlock(_disp_x_var, _disp_z_var);
//    DenseMatrix<Number> & A10 = _assembly.jacobianBlock(_disp_y_var, _disp_x_var);
//    DenseMatrix<Number> & A11 = _assembly.jacobianBlock(_disp_y_var, _disp_y_var);
//    DenseMatrix<Number> & A12 = _assembly.jacobianBlock(_disp_y_var, _disp_z_var);
//    DenseMatrix<Number> & A20 = _assembly.jacobianBlock(_disp_z_var, _disp_x_var);
//    DenseMatrix<Number> & A21 = _assembly.jacobianBlock(_disp_z_var, _disp_y_var);
//    DenseMatrix<Number> & A22 = _assembly.jacobianBlock(_disp_z_var, _disp_z_var);
//
//    DenseMatrix<Number> & B03 = _assembly.jacobianBlock(_disp_x_var, _pres_var  );
//    DenseMatrix<Number> & B13 = _assembly.jacobianBlock(_disp_y_var, _pres_var  );
//    DenseMatrix<Number> & B23 = _assembly.jacobianBlock(_disp_z_var, _pres_var  );
//    DenseMatrix<Number> & B30 = _assembly.jacobianBlock(_pres_var  , _disp_x_var);
//    DenseMatrix<Number> & B31 = _assembly.jacobianBlock(_pres_var  , _disp_y_var);
//    DenseMatrix<Number> & B32 = _assembly.jacobianBlock(_pres_var  , _disp_z_var);
//
//
//    DenseMatrix<Number> _local_A[3][3];
//    DenseMatrix<Number> _local_BT[3];
//    DenseMatrix<Number> _local_B [3];
//
//    // set dimensions of local Jacobians and set them to zero
//    for (int i=0; i<3; ++i)
//        for (int j=0; j<3; ++j)
//        {
//            _local_A[i][j].resize(A00.m(),A00.n());
//            _local_A[i][j].zero();
//        }
//
//
//    if(_pres_var!= NULL){
//
//        for (int j=0; j<3; ++j)
//        {
//            _local_BT[j].resize(B03.m(),B03.n());
//            _local_BT[j].zero();
//        }
//
//
//        for (int j=0; j<3; ++j)
//        {
//            _local_B[j].resize(B30.m(),B30.n());
//            _local_B[j].zero();
//        }
//    }
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
//                if(_pres_var!= NULL)_stress_lin[dim][j][qp]+=_p[qp]*_J[qp]* ( _invFtr[qp].contract(H)*_invFtr[qp]-_invFtr[qp]*H.transpose()*_invFtr[qp]  );
//                _stress_lin[dim][j][qp] = 1./4 *_stress_lin[dim][j][qp];
//            }
//        }
//    }
//
//    for (int idim=0; idim<3; ++idim)
//        for (int jdim=idim; jdim<3; ++jdim)
//            for (_i = 0; _i < _test.size(); _i++)
//                for (_j = 0; _j < _phi.size(); _j++)
//                    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
//                        _local_A[idim][jdim](_i,_j)+=_JxW[_qp] * _coord[_qp]*
//                        (_stress_lin[jdim][_j][_qp](idim,0)*_grad_test[_i][_qp](0)+
//                         _stress_lin[jdim][_j][_qp](idim,1)*_grad_test[_i][_qp](1)+ // forse jdim?
//                         _stress_lin[jdim][_j][_qp](idim,2)*_grad_test[_i][_qp](2)
//                         );
//
//
//    if(_pres_var!= NULL){
//        for (int idim=0; idim<3; ++idim)
//            for (_i=0; _i<_test.size(); ++_i)
//                for (_j=0; _j<_test_pres->size(); ++_j)
//                    for (_qp=0; _qp<_qrule->n_points(); ++_qp)
//                        _local_BT[idim](_i,_j)+=1./4*_JxW[_qp] * _coord[_qp]*((_test_pres[0][_j][_qp]))*_J[_qp]*
//                        (_invFtr[_qp](idim,0)*_grad_test[_i][_qp](0)+_invFtr[_qp](idim,1)*_grad_test[_i][_qp](1)+_invFtr[_qp](idim,2)*_grad_test[_i][_qp](2));
//    }
//
//
//    _local_A[0][1].get_transpose(_local_A[1][0]);
//    _local_A[0][2].get_transpose(_local_A[2][0]);
//    _local_A[1][2].get_transpose(_local_A[2][1]);
//
//    if(_pres_var!= NULL){
//        _local_BT[0].get_transpose(_local_B[0]);
//        _local_BT[1].get_transpose(_local_B[1]);
//        _local_BT[2].get_transpose(_local_B[2]);
//    }
//
//    A00 += _local_A[0][0];
//    A01 += _local_A[0][1];
//    A02 += _local_A[0][2];
//
//    A10 += _local_A[1][0];
//    A11 += _local_A[1][1];
//    A12 += _local_A[1][2];
//
//    A20 += _local_A[2][0];
//    A21 += _local_A[2][1];
//    A22 += _local_A[2][2];
//
//    if(_pres_var!= NULL){
//
//        B03 += _local_BT[0];
//        B13 += _local_BT[1];
//        B23 += _local_BT[2];
//
//        B30 += _local_B[0];
//        B31 += _local_B[1];
//        B32 += _local_B[2];
//
//    }
//
//
//
//    if(_has_diag_save_in)
//    {
//        mooseError("Error: diag in already saved in");
//    }
//
//}
//
//
//void NewmarkStressDivergenceFast::computeResidual()
//{
//    if (_dim==2)
//        computeResidual2D();
//    else
//    {
//        if (_dim==3)
//            computeResidual3D();
//        else
//            mooseError( "You should not be here");
//
//    }
//}
//
//void
//NewmarkStressDivergenceFast::computeResidual2D()
//{
//    // References to the residual components
//    DenseVector<Number> & re_x = _assembly.residualBlock(_disp_x_var);
//    DenseVector<Number> & re_y = _assembly.residualBlock(_disp_y_var);
//    DenseVector<Number> & re_p = _assembly.residualBlock(_pres_var);
//
//    // Local vectors where the local contribution is assembled
//    DenseVector<Number> _local_re[2];
//    for (int i=0; i<2; ++i)
//    {
//        _local_re[i].resize(re_x.size());
//        _local_re[i].zero();
//    }
//
//    DenseVector<Number> _local_re_p;
//    _local_re_p.resize(re_p.size());
//
//    // Assembly of residual for mechanics
//    for (int d=0; d<2; ++d)
//        for (_i = 0; _i < _test.size(); ++_i)
//            for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
//            {
//                if (_t_step==1){
//
//                    RealTensorValue P= 1./4 * ( _P[_qp] +_P_old[_qp]);
//
//                    if(_pres_var!= NULL) P +=  1./4 * ( _p[_qp]*_J[_qp] *_invFtr[_qp] + _p_old[_qp]*_J_old[_qp]*_invFtr_old[_qp]);
//
//                    _local_re[d](_i) += _JxW[_qp] * _coord[_qp] *
//                    (
//                     P(d , 0) * _grad_test[_i][_qp] (0)+
//                     P(d , 1) * _grad_test[_i][_qp] (1)
//                     );
//                }
//                if (_t_step>1){
//
//                    RealTensorValue P= 1./4 * ( _P[_qp] ) +
//                                       1./2 * (_P_old[_qp])+
//                                       1./4 * (_P_older[_qp]) ;
//
//                    if(_pres_var!= NULL) {P += 1./4 * (_p[_qp]*_J[_qp] *_invFtr[_qp]) +
//                                                1./2 * (_p_old[_qp]*_J_old[_qp]*_invFtr_old[_qp])+
//                                                1./4 * (_p_older[_qp] * _J_older[_qp] * _invFtr_older[_qp]) ;}
//
//                    _local_re[d](_i) += _JxW[_qp] * _coord[_qp] *
//                                                                (
//                                                                P(d , 0) * _grad_test[_i][_qp] (0)+
//                                                                P(d , 1) * _grad_test[_i][_qp] (1)
//                                                                );
//                }
//            }
//
//    if(_pres_var!= NULL) {
//        for (_i =0; _i<_test_pres->size(); ++_i)
//            for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
//                _local_re_p(_i)+=1./4 * _JxW[_qp] * _coord[_qp] *(_J[_qp]-_J_old[_qp])* (_test_pres[0][_i][_qp]);
//    }
//
//    // sum the local contributions to the global residual
//    re_x += _local_re[0];
//    re_y += _local_re[1];
//    if(_pres_var!= NULL) re_p += _local_re_p;
//
//}
//
//void
//NewmarkStressDivergenceFast::computeResidual3D()
//{
//    // References to the residual components
//    DenseVector<Number> & re_x = _assembly.residualBlock(_disp_x_var);
//    DenseVector<Number> & re_y = _assembly.residualBlock(_disp_y_var);
//    DenseVector<Number> & re_z = _assembly.residualBlock(_disp_z_var);
//    DenseVector<Number> & re_p = _assembly.residualBlock(_pres_var);
//
//    // Local vectors where the local contribution is assembled
//    DenseVector<Number> _local_re[3];
//    for (int i=0; i<3; ++i)
//    {
//        _local_re[i].resize(re_x.size());
//        _local_re[i].zero();
//    }
//
//    DenseVector<Number> _local_re_p;
//    _local_re_p.resize(re_p.size());
//
//    // Assembly of residual for mechanics
//    for (int d=0; d<3; ++d)
//        for (_i = 0; _i < _test.size(); ++_i)
//            for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
//            {
//                if (_t_step==1){
//
//                    RealTensorValue P= 1./4 * ( _P[_qp] +_P_old[_qp]);
//
//                    if(_pres_var!= NULL) P +=  1./4 * ( _p[_qp]*_J[_qp] *_invFtr[_qp] + _p_old[_qp]*_J_old[_qp]*_invFtr_old[_qp]);
//
//                    _local_re[d](_i) += _JxW[_qp] * _coord[_qp] *
//                    (
//                     P(d , 0) * _grad_test[_i][_qp] (0)+
//                     P(d , 1) * _grad_test[_i][_qp] (1)+
//                     P(d , 2) * _grad_test[_i][_qp] (2));
//                }
//                if (_t_step>1){
//
//                    RealTensorValue P= 1./4 * ( _P[_qp] ) +
//                    1./2 * (_P_old[_qp])+
//                    1./4 * (_P_older[_qp]) ;
//
//                    if(_pres_var!= NULL) {P += 1./4 * (_p[_qp]*_J[_qp] *_invFtr[_qp]) +
//                        1./2 * (_p_old[_qp]*_J_old[_qp]*_invFtr_old[_qp])+
//                        1./4 * (_p_older[_qp] * _J_older[_qp] * _invFtr_older[_qp]) ;}
//
//                    _local_re[d](_i) += _JxW[_qp] * _coord[_qp] *
//                    (
//                     P(d , 0) * _grad_test[_i][_qp] (0)+
//                     P(d , 1) * _grad_test[_i][_qp] (1)+
//                     P(d , 2) * _grad_test[_i][_qp] (2));
//                }
//            }
//
//
//
//    if(_pres_var!= NULL) {
//        for (_i =0; _i<_test_pres->size(); ++_i)
//            for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
//                _local_re_p(_i)+=_JxW[_qp] * _coord[_qp] *(_J[_qp]-1)* (_test_pres[0][_i][_qp]);
//    }
//
//    // sum the local contributions to the global residual
//    re_x += _local_re[0];
//    re_y += _local_re[1];
//    re_z += _local_re[2];
//    if(_pres_var!= NULL) re_p += _local_re_p;
//
//}
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
/*                Prepared by Maria Nestola                     */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/****************************************************************/



#include "NewmarkStressDivergenceFast.h"
#include "Material.h"
#include "Assembly.h"
#include "MooseMesh.h"
#include "SystemBase.h"
#include "MooseVariable.h"
#include "libmesh/quadrature.h"

registerMooseObject("MECHApp", NewmarkStressDivergenceFast);

template <>
InputParameters validParams<NewmarkStressDivergenceFast>() {
    
    // inherit the parameters of the Kernels:
    InputParameters params = validParams<Kernel>();
    params.addRequiredCoupledVar("disp_y","disp_y");
    params.addCoupledVar("disp_z","disp_z");
    //params.addCoupledVar("pressure","pressure");
    params.addParam<NonlinearVariableName>("pressure", "pressure");
    params.addParam<bool>("linear",false, "linear");
    
    
    return params;
}


NewmarkStressDivergenceFast::NewmarkStressDivergenceFast(
                                           InputParameters const & parameters) :

// inherit the parameters of the Kernels:
Kernel(parameters),
// We are doing elasticity from R^n to R^n,
// hence the number of components is the same as the mesh dimension:
// For example for _mesh.dimension=3 we have
//_component = 0 for disp_x, _component = 1 for disp_y, _component = 2 for
// disp_z
_dim(_mesh.dimension()),
_disp_x_var(_var.number()),
_disp_y_var(coupled("disp_y")),
_disp_z_var(_mesh.dimension() == 3 ? coupled("disp_z") : 100000),
_has_pressure(parameters.isParamValid("pressure") ? 1 : 0),
_is_linear(getParam<bool>("linear")),
// inherit material properties:
_p(getMaterialProperty<Real>("pressureFromMaterial")),
_J(getMaterialProperty<Real>("deformationDeterminant")),
_invFtr(getMaterialProperty<RealTensorValue>("invFtr")),
_P(getMaterialProperty<RealTensorValue>("firstPiolaKirchhofTensor")),
_P_old(getMaterialPropertyOld<RealTensorValue>("firstPiolaKirchhofTensor")),
_P_older(getMaterialPropertyOlder<RealTensorValue>("firstPiolaKirchhofTensor")),
_materialPointer(getMaterialProperty<ElasticMaterial *>("materialPointer")),
_gradDisp(getMaterialProperty<RealTensorValue>("displacementGradient")),
_variable_pres(NULL),
_test_pres(NULL)

{
    // we suppose that the components of displacements are the first three.
    
    if(_has_pressure)
    {
        _variable_pres = &_subproblem.getStandardVariable(_tid, getParam<NonlinearVariableName>("pressure"));
        _test_pres = &(_variable_pres[0].phi());
        _pres_var=_variable_pres[0].number() ;
    }
    
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
    
    _residual=new DenseVector<Number> * [_dim];
    
    _local_re=new DenseVector<Number> [_dim];
    
    _disp_var=new unsigned [_dim];
    
    if (_dim==2)
    {
        _disp_var[0]=_disp_x_var;
        _disp_var[1]=_disp_y_var;
    }
    
    if (_dim==3)
    {
        _disp_var[0]=_disp_x_var;
        _disp_var[1]=_disp_y_var;
        _disp_var[2]=_disp_z_var;
    }
    
    _jacobian=new DenseMatrix<Number> ** [_dim];
    for (int d=0; d<_dim; ++d)
        _jacobian[d]=new DenseMatrix<Number> * [_dim];
    
    _local_A=new DenseMatrix<Number> * [_dim];
    for (int d=0; d<_dim; ++d)
        _local_A[d]=new DenseMatrix<Number> [_dim];
    
    if (_has_pressure)
    {
        _B        =new DenseMatrix<Number> * [_dim];
        _local_B  =new DenseMatrix<Number>   [_dim];
        _BT       =new DenseMatrix<Number> * [_dim];
        _local_BT =new DenseMatrix<Number>   [_dim];
    }
    
}

NewmarkStressDivergenceFast::~NewmarkStressDivergenceFast() {
    for (int d=0; d<_dim; ++d) {
        for (int j = 0; j < 64 /*_phi.size()*/; ++j) {
            delete[] _stress_lin[d][j];
        }
        delete[] _stress_lin[d];
        delete[] _jacobian[d];
        delete[] _local_A[d];
    }
    delete[] _stress_lin;
    delete[] _residual;
    delete[] _local_re;
    delete[] _disp_var;
    delete[] _jacobian;
    delete[] _local_A;
    
    if (_has_pressure)
    {
        delete[] _B;
        delete[] _local_B;
        delete[] _BT;
        delete[] _local_BT;
    }
}

void
NewmarkStressDivergenceFast::computeOffDiagJacobian(unsigned int jvar)
{
    if (jvar==_disp_x_var)
        computeJacobian();
}




void
NewmarkStressDivergenceFast::computeJacobian()
{
    

    for (int idim=0; idim<_dim; ++idim)
        for (int jdim=0; jdim<_dim; ++jdim)
        {
            _jacobian[idim][jdim]=&_assembly.jacobianBlock(_disp_var[idim], _disp_var[jdim]);
            _local_A[idim][jdim].resize(_test.size(),_test.size());
            _local_A[idim][jdim].zero();
        }
    
    
    for (int jdim=0; jdim<_dim; ++jdim)
    {
        H = RealTensorValue(zeros, zeros, zeros);
        for (_j=0; _j<_phi.size(); ++_j)
        {
            for (int _qp=0; _qp<_qrule->n_points(); ++_qp)
            {
                for (int k=0 ; k<_dim; ++k)
                    H(jdim,k)=_grad_phi[_j][_qp](k);
                (*(_materialPointer[_qp])).evaluate_stress_lin(_qp, H, _stress_lin[jdim][_j][_qp]);
                if(_has_pressure) _stress_lin[jdim][_j][_qp]+=_p[_qp]*_J[_qp]* ( _invFtr[_qp].contract(H)*_invFtr[_qp]-_invFtr[_qp]*H.transpose()*_invFtr[_qp]  );
            }
        }
    }
    
    for (int idim=0; idim<_dim; ++idim)
        for (int jdim=idim; jdim<_dim; ++jdim)
            for (_i = 0; _i < _test.size(); _i++)
                for (_j = 0; _j < _phi.size(); _j++)
                    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
                        for (int k=0; k<_dim; ++k)
                            _local_A[idim][jdim](_i,_j)+=_JxW[_qp] * _coord[_qp]*_stress_lin[jdim][_j][_qp](idim,k)*_grad_test[_i][_qp](k);
    
    
    for (int idim=1; idim<_dim; ++idim)
        for (int jdim=0; jdim<idim; ++jdim)
            _local_A[jdim][idim].get_transpose(_local_A[idim][jdim]);
    
    
    for (int idim=0; idim<_dim; ++idim)
        for (int jdim=0; jdim<_dim; ++jdim)
        {
            _jacobian[idim][jdim][0]+=_local_A[idim][jdim];
        }
    
    
    
    if (_has_pressure)
    {
        for (int idim=0; idim<_dim; ++idim)
        {
            _BT[idim]=&_assembly.jacobianBlock(_disp_var[idim], _pres_var  );
            _local_BT[idim].resize(_test.size(),_test_pres[0].size());
            _local_BT[idim].zero();

            _B[idim]=&_assembly.jacobianBlock(_pres_var  , _disp_var[idim]);
            _local_B[idim].resize(_test_pres[0].size(),_test.size());
            _local_B[idim].zero();
        }


        for (int idim=0; idim<_dim; ++idim){
            for (_i=0; _i<_test.size(); ++_i){
                for (_j=0; _j<_test_pres[0].size(); ++_j){
                    for (_qp=0; _qp<_qrule->n_points(); ++_qp){
                        for (int k=0; k<_dim; ++k){
                            if (_is_linear) _local_BT[idim](_i,_j)+=_JxW[_qp] * _coord[_qp]*(_test_pres[0][_j][_qp])*(1 + _gradDisp[_qp].tr())*_invFtr[_qp](idim,k)*_grad_test[_i][_qp](k);
                            else _local_BT[idim](_i,_j)+=_JxW[_qp] * _coord[_qp]*(_test_pres[0][_j][_qp])*_J[_qp]*_invFtr[_qp](idim,k)*_grad_test[_i][_qp](k);
                        }
                    }
                }
            }
        }


        for (int d=0; d<_dim; ++d)
            _local_BT[d].get_transpose(_local_B[d]);

        for (int d=0; d<_dim; ++d)
        {
            _BT[d][0] += _local_BT[d];
            _B[d][0] += _local_B[d];
        }




    }
    
    
    
    if(_has_diag_save_in)
    {
        mooseError("Error: diag in already saved in");
    }
    
    
    
}


void NewmarkStressDivergenceFast::computeResidual()
{
    //    if (_current_elem->id()==1)
    //    {
    //        for (_qp=0; _qp < _qrule->n_points(); ++_qp)
    //             std::cout<<"RES "<<_q_point[_qp]<<" p="<<_p[_qp]<<std::endl;
    //
    //             std::cout<<std::endl<<std::endl;
    //    }
    for (int d=0; d<_dim; ++d)
    {
        _residual[d] = &_assembly.residualBlock(_disp_var[d]);
        _local_re[d].resize(_test.size());
        _local_re[d].zero();
    }
    
    // Assembly of residual for mechanics
    for (int d=0; d<_dim; ++d)
        for (_i = 0; _i < _test.size(); ++_i)
            for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
            {
                RealTensorValue  P = _P[_qp];
                
                if (_t_step==1){
                              P += _P_old[_qp];
                }
                else{
                              P += 2 * _P_old[_qp] +  _P_older[_qp];
                }
                
                for (int k=0; k<_dim; ++k)
                    _local_re[d](_i) += _JxW[_qp] * _coord[_qp] * P(d , k) * _grad_test[_i][_qp] (k);
            }
    
    
    
    for (int d=0 ; d<_dim; ++d)
    {
        _residual[d][0] += _local_re[d];
    }
    
    
    
    //
    //    if (_has_save_in)
    //    {
    //        Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    //        for (const auto & var : _save_in)
    //            var->sys().solution().add_vector(_residual[0][0], var->dofIndices());
    //    }
    
    
    //  DenseVector<Number> & re_p = _assembly.residualBlock(_pres_var);
    
    
    if(_has_pressure)
    {
        // DenseVector<Number> & re_p = _assembly.residualBlock(_pres_var);//, Moose::KT_TIME);
        
        DenseVector<Number> & re_p = _assembly.residualBlock(_pres_var);
        DenseVector<Number> _local_re_p;
        _local_re_p.resize(_test_pres[0].size());
        
        for (_i =0; _i<_test_pres->size(); ++_i){
            for (_qp = 0; _qp < _qrule->n_points(); ++_qp){
                if (_is_linear) _local_re_p(_i)+=_JxW[_qp] * _coord[_qp] * _gradDisp[_qp].tr() * (_test_pres[0][_i][_qp]);
                else _local_re_p(_i)+=_JxW[_qp] * _coord[_qp] *(_J[_qp]-1.0)* (_test_pres[0][_i][_qp]);
            }
        }
        
        re_p += _local_re_p;
        
        
    }
    
    
    if (_has_save_in){
        
        for (unsigned int i=0; i<_save_in_strings.size(); i++)
        {
            Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
            MooseVariable * var = &_subproblem.getStandardVariable(_tid, _save_in_strings[i]);
            var->sys().solution().add_vector( _local_re[i],var->dofIndices());
            //std::cout << _local_re[0].size() << std::endl;
            //std::cout<<"_local_re[0] = "<<_local_re[0] <<std::endl;
            
        }
    }
    
    
    
    
}
