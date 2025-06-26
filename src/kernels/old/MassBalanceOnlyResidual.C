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
/* Kernel to ensure the mass balance. To be used only coupled   */
/* with IncompressibleStressDivergenceLibmeshTensor. Only       */
/* residuals computation.                                       */
/****************************************************************/


#include "MassBalanceOnlyResidual.h"
#include "MooseMesh.h"

#include "Assembly.h"
#include "MooseVariable.h"
#include "libmesh/quadrature.h"

registerMooseObject("MECHApp", MassBalanceOnlyResidual);

template <>
InputParameters validParams<MassBalanceOnlyResidual>() {

 // inherit the parameters of the Kernels:
   InputParameters params = validParams<Kernel>();

   // specify for which component to solve:
  params.addRequiredParam<unsigned>("component", "component");
  return params;
}


MassBalanceOnlyResidual::MassBalanceOnlyResidual(
                                         InputParameters const & parameters) :

  // inherit the parameters of the Kernels:
    Kernel(parameters),

  // specify for which component to solve:
    _component(getParam<unsigned>("component")),
 
  // We are doing elasticity from R^n to R^n,
  // hence the number of components is the same as the mesh dimension:
  // For example for _mesh.dimension=3 we have
  //_component = 0 for disp_x, _component = 1 for disp_y, _component = 2 for
  // disp_z
    _numComp(_mesh.dimension()),
  
  // inherit some material properties:
    _J(getMaterialProperty<Real>("deformationDeterminant"))
{}

Real MassBalanceOnlyResidual::computeQpResidual() {

    return (_J[_qp]-1.0)*_test[_i][_qp];
}


void
MassBalanceOnlyResidual::computeJacobian()
{
    DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
    
    if(_has_diag_save_in)
    { mooseError("Error: diag in already saved in");
    }
}

void
MassBalanceOnlyResidual::computeOffDiagJacobian(unsigned int jvar)
{
    
    DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), jvar);
    
    if(_has_diag_save_in)
    { mooseError("Error: diag in already saved in");
    }
    
}
