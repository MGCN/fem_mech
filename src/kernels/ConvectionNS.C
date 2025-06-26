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



#include "ConvectionNS.h"

#include "Material.h"
#include "Assembly.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "SystemBase.h"

#include "libmesh/quadrature.h"


registerMooseObject("MECHApp", ConvectionNS);

template <>
InputParameters validParams<ConvectionNS>() {

  // inherit the parameters of the Kernels:
  InputParameters params = validParams<Kernel>();
  params.addRequiredCoupledVar("vel_y","vel_y");
  params.addCoupledVar("vel_z","vel_z");
  
  return params;
}


ConvectionNS::ConvectionNS(InputParameters const & parameters) :

    // inherit the parameters of the Kernels:
    Kernel(parameters),
    // We are doing elasticity from R^n to R^n,
    // hence the number of components is the same as the mesh dimension:
    // For example for _mesh.dimension=3 we have
    //_component = 0 for disp_x, _component = 1 for disp_y, _component = 2 for
    // disp_z
    _dim(_mesh.dimension()),
    _vel_x_var(_var.number()),
    _vel_y_var(coupled("vel_y")),
    _vel_z_var(_mesh.dimension() == 3 ? coupled("vel_z") : 100000),
    _grad_vel(getMaterialProperty<RealTensorValue>("displacementGradient")),
    _vel_y(coupledValue("vel_y")),
    _vel_z(_mesh.dimension() == 3 ? coupledValue("vel_z") : _zero)
{

  if (_dim == 1)
  {
     mooseError( "Stress divergence fast does not work with dim=1");
  }
    
    _stress_lin = new RealVectorValue **[_dim] /*[_phi.size()]*/;
    
    for (int d=0; d<_dim; ++d)
        _stress_lin[d]=new RealVectorValue *[64];
    
    
    for (int d=0; d<_dim; ++d)
        for (int j = 0; j < 64 /*_phi.size()*/; ++j)
            _stress_lin[d][j] = new RealVectorValue [64] /*[_qrule->n_points()]*/;
    
    _residual=new DenseVector<Number> * [_dim];
    
    _local_re=new DenseVector<Number> [_dim];
    
    _vel_var=new unsigned [_dim];
    
    if (_dim==2)
    {
        _vel_var[0]=_vel_x_var;
        _vel_var[1]=_vel_y_var;
    }

    if (_dim==3)
    {
        _vel_var[0]=_vel_x_var;
        _vel_var[1]=_vel_y_var;
        _vel_var[2]=_vel_z_var;
    }

    _jacobian=new DenseMatrix<Number> ** [_dim];
    for (int d=0; d<_dim; ++d)
        _jacobian[d]=new DenseMatrix<Number> * [_dim];
    
    _local_A=new DenseMatrix<Number> * [_dim];
    for (int d=0; d<_dim; ++d)
        _local_A[d]=new DenseMatrix<Number> [_dim];
    
    _vel=new RealVectorValue [64];
    _non_linear_term=new RealVectorValue [64];

}

void
ConvectionNS::computeOffDiagJacobian(unsigned int jvar)
{
    if (jvar==_vel_x_var)
        computeJacobian();
}




void
ConvectionNS::computeJacobian()
{
    
    for (int idim=0; idim<_dim; ++idim)
        for (int jdim=0; jdim<_dim; ++jdim)
        {
            _jacobian[idim][jdim]=&_assembly.jacobianBlock(_vel_var[idim], _vel_var[jdim]);
            _local_A[idim][jdim].resize(_test.size(),_test.size());
            _local_A[idim][jdim].zero();
        }
    
    
    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
    {
        _vel[_qp](0)=_u[_qp];
        _vel[_qp](1)=_vel_y[_qp];
        _vel[_qp](2)=_vel_z[_qp];
    }

    
    for (int jdim=0; jdim<_dim; ++jdim)
    {
        _H = RealTensorValue(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
        _h = RealVectorValue(0.0,0.0,0.0);
        for (_j=0; _j<_phi.size(); ++_j)
        {
            for (int _qp=0; _qp<_qrule->n_points(); ++_qp)
            {
                for (int k=0 ; k<_dim; ++k)
                    _H(jdim,k)=_grad_phi[_j][_qp](k);
                
                _h(jdim)=_phi[_j][_qp];
                
                RealTensorValue const & U=_grad_vel[_qp];
                RealVectorValue const & u=_vel[_qp];
                _stress_lin[jdim][_j][_qp]=_H*u+U*_h;

            }
        }
    }
    
    for (int idim=0; idim<_dim; ++idim)
        for (int jdim=idim; jdim<_dim; ++jdim)
            for (_i = 0; _i < _test.size(); _i++)
                for (_j = 0; _j < _phi.size(); _j++)
                    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
                            _local_A[idim][jdim](_i,_j)+=_JxW[_qp] * _coord[_qp]*_stress_lin[jdim][_j][_qp](idim)*_test[_i][_qp];

    
    for (int idim=1; idim<_dim; ++idim)
        for (int jdim=0; jdim<idim; ++jdim)
            _local_A[jdim][idim].get_transpose(_local_A[idim][jdim]);

    
    for (int idim=0; idim<_dim; ++idim)
        for (int jdim=0; jdim<_dim; ++jdim)
        {
            _jacobian[idim][jdim][0]+=_local_A[idim][jdim];
        }

    
    if(_has_diag_save_in)
    {
        mooseError("Error: diag in already saved in");
    }


    
}


void ConvectionNS::computeResidual()
{
    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
    {
        _vel[_qp](0)=_u[_qp];
        _vel[_qp](1)=_vel_y[_qp];
        _vel[_qp](2)=_vel_z[_qp];
        
        RealTensorValue const & U=_grad_vel[_qp];
        RealVectorValue const & u=_vel[_qp];
        _non_linear_term[_qp]=U*u;
    }
    
    for (int d=0; d<_dim; ++d)
    {
        _residual[d] = &_assembly.residualBlock(_vel_var[d]);
        _local_re[d].resize(_test.size());
        _local_re[d].zero();
    }

    // Assembly of residual for mechanics
    for (int d=0; d<_dim; ++d)
        for (_i = 0; _i < _test.size(); ++_i)
            for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
            {
                RealVectorValue const & NL = _non_linear_term[_qp];
                
                    _local_re[d](_i) += _JxW[_qp] * _coord[_qp] * NL(d) * _test[_i][_qp];
            }
    
    

    for (int d=0 ; d<_dim; ++d)
    {
        _residual[d][0] += _local_re[d];
    }
    
    
    for (unsigned int i=0; i<_save_in_strings.size(); i++)
    {

      MooseVariable * var = &_subproblem.getStandardVariable(_tid, _save_in_strings[i]);
      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
      var->sys().solution().add_vector(_residual[i][0],var->dofIndices());

    }


    
}
