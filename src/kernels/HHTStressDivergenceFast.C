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



#include "HHTStressDivergenceFast.h"

#include "Material.h"
#include "Assembly.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "SystemBase.h"

#include "libmesh/quadrature.h"
registerMooseObject("MECHApp", HHTStressDivergenceFast);

template <>
InputParameters validParams<HHTStressDivergenceFast>() {
    
    InputParameters params = validParams<Kernel>();
    params.addRequiredCoupledVar("disp_y","disp_y");
    params.addCoupledVar("disp_z","disp_z");
    params.addCoupledVar("pressure","pressure");
    params.addParam<NonlinearVariableName>("variable_pres", "Pressure variable");
    params.addParam<Real>("alpha",false,
                          "HHT time inetgartion");
    return params;
}


HHTStressDivergenceFast::HHTStressDivergenceFast(InputParameters const & parameters) :

// the pressure variable in coupled and non coupled copy:
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
_pres_var(parameters.isParamValid("pressure") ? coupled("pressure") : 100000 ),
_has_pressure(parameters.isParamValid("pressure") ? 1 : 0),
_p(parameters.isParamValid("pressure") ? coupledValue("pressure"):_zero),
// inherit material properties:
_J(getMaterialProperty<Real>("deformationDeterminant")),
_invFtr(getMaterialProperty<RealTensorValue>("invFtr")),
_P(getMaterialProperty<RealTensorValue>("firstPiolaKirchhofTensor")),
_P_old(getMaterialPropertyOld<RealTensorValue>("firstPiolaKirchhofTensor")),
_materialPointer(getMaterialProperty<ElasticMaterial *>("materialPointer")),
_alpha_f(getParam<Real>("alpha")),
_variable_pres(NULL),
_test_pres(NULL)


{
    
    if(_has_pressure){
        _variable_pres = &_subproblem.getStandardVariable(_tid, "variable_pres");
        _test_pres = &_variable_pres->phi();
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

void
HHTStressDivergenceFast::computeOffDiagJacobian(unsigned int jvar)
{
    if (jvar==_disp_x_var)
        computeJacobian();
}


void
HHTStressDivergenceFast::computeJacobian()
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
                if(_has_pressure)
                    _stress_lin[jdim][_j][_qp]+=_p[_qp]*_J[_qp]* ( _invFtr[_qp].contract(H)*_invFtr[_qp]-_invFtr[_qp]*H.transpose()*_invFtr[_qp]  );
                    _stress_lin[jdim][_j][_qp]=( 1.0 - _alpha_f ) *_stress_lin[jdim][_j][_qp];
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
        
        
        for (int idim=0; idim<_dim; ++idim)
            for (_i=0; _i<_test.size(); ++_i)
                for (_j=0; _j<_test_pres[0].size(); ++_j)
                    for (_qp=0; _qp<_qrule->n_points(); ++_qp)
                        for (int k=0; k<_dim; ++k)
                            _local_BT[idim](_i,_j)+=_JxW[_qp] * _coord[_qp]*(_test_pres[0][_j][_qp])*_J[_qp]*_invFtr[_qp](idim,k)*_grad_test[_i][_qp](k);
        
        
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



void HHTStressDivergenceFast::computeResidual()
{
    for (int d=0; d<_dim; ++d)
    {
        _residual[d] = &_assembly.residualBlock(_disp_var[d]);
        _local_re[d].resize(_test.size());
        _local_re[d].zero();
    }

    RealVectorValue _reaction(0.0,0.0,0.0);
    
    // Assembly of residual for mechanics
    for (int d=0; d<_dim; ++d)
        for (_i = 0; _i < _test.size(); ++_i)
            for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
            {
                     RealTensorValue const & P = (1.0 - _alpha_f) * _P[_qp] + _alpha_f * _P_old[_qp];
                    
                    for (int k=0; k<_dim; ++k){
                        _local_re[d](_i) += _JxW[_qp] * _coord[_qp] * P(d , k) * _grad_test[_i][_qp] (k);

                     }
                
            }
    
    for (int d=0 ; d<_dim; ++d)
        _residual[d][0] += _local_re[d];
    
    
    if(_has_pressure)
    {
        DenseVector<Number> & re_p = _assembly.residualBlock(_pres_var);
        
        DenseVector<Number> _local_re_p;
        _local_re_p.resize(_test_pres[0].size());
        
        for (_i =0; _i<_test_pres->size(); ++_i)
            for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
                _local_re_p(_i)+= _JxW[_qp] * _coord[_qp] *(_J[_qp]-1.0)* (_test_pres[0][_i][_qp]);
        
        re_p += _local_re_p;
        
        
    }


    if (_has_save_in){
        
        for (unsigned int i=0; i<_save_in_strings.size(); i++)
        {
            Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
            MooseVariable * var = &_subproblem.getStandardVariable(_tid, _save_in_strings[i]);
            var->sys().solution().add_vector( _local_re[i],var->dofIndices());
            //std::cout << _local_re[0].size() << std::endl;
            // std::cout<<"_local_re[0] = "<<_local_re[0] <<std::endl;
            
        }
    }
    
    
}


