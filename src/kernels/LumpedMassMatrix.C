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
/****************************************************************/


#include "LumpedMassMatrix.h"
#include "MooseVariable.h"
#include "Assembly.h"
#include "MooseMesh.h"
#include "libmesh/quadrature.h"

registerMooseObject("MECHApp", LumpedMassMatrix);
template<>
InputParameters validParams<LumpedMassMatrix>()
{
  InputParameters params = validParams<Kernel>();
  return params;
}

LumpedMassMatrix::LumpedMassMatrix(const InputParameters & parameters) :
    Kernel(parameters),
    _u_nodal(_var.dofValues())
{}

Real
LumpedMassMatrix::computeQpResidual()
{

  return 0.0;
}

Real
LumpedMassMatrix::computeQpJacobian()
{
  return 0.0;
}


void
LumpedMassMatrix::computeResidual()
{
       
    DenseVector<Number> & re = _assembly.residualBlock(_var.number());
    
    DenseVector<Number> _local_re;
    
    _local_re.resize(re.size());
    
    _local_re.zero();

    precalculateResidual();
          
    for (_i = 0; _i < _test.size(); _i++)
        for (_qp = 0; _qp < _qrule->n_points(); _qp++)
               _local_re(_i) += _JxW[_qp] * _coord[_qp] * _u_nodal[_i] * _test[_i][_qp];

        re += _local_re;
    
    
    
}

void
LumpedMassMatrix::computeJacobian()
{
    

    DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
    _local_ke.resize(ke.m(), ke.n());
    _local_ke.zero();

    
    precalculateJacobian();
    for (_i = 0; _i < _test.size(); _i++)
        for (_j = 0; _j < _phi.size(); _j++)
            for (_qp = 0; _qp < _qrule->n_points(); _qp++){
                _local_ke(_i, _j) += _JxW[_qp] * _coord[_qp] * _test[_i][_qp]*_phi[_j][_qp];
                //std::cout<<"qpoints = "<<_qrule->n_points()<<std::endl;
            }
    
    //std::cout<<_local_ke <<std::endl<<std::endl;
    
    for (int i=0; i<ke.m(); ++i)
        for (int j=0; j<ke.n(); ++j)
            if (i!=j)
            {
                _local_ke(i,i)+=_local_ke(i,j);
                _local_ke(i,j)=0.0;
            }
//     std::cout<<_local_ke<<std::endl<<std::endl;
    
    
    ke += _local_ke;
    
    
}
