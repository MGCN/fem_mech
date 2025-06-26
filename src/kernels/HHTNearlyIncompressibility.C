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
/*                                                              */
/*         Kernel for the Nearly Incompressibility              */
/****************************************************************/

#include "HHTNearlyIncompressibility.h"
#include "libmesh/quadrature.h"
#include "MooseMesh.h"
#include "MooseVariable.h"

registerMooseObject("MECHApp", HHTNearlyIncompressibility);

template <>
InputParameters validParams<HHTNearlyIncompressibility>() {
InputParameters params = validParams<Kernel>();
params.addRequiredCoupledVar("displacements", "The displacements appropriate for the simulation geometry and coordinate system");
params.addRequiredParam<Real>("kappa","bulk_modulus");
params.addRequiredParam<Real>("alpha","HHT parameter");
return params;
}


HHTNearlyIncompressibility::HHTNearlyIncompressibility(InputParameters const & parameters) :
    Kernel(parameters),
    _numComp(_mesh.dimension()),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_numComp),
    _J(getMaterialProperty<Real>("deformationDeterminant")),
    _J_old(getMaterialPropertyOld<Real>("deformationDeterminant")),
    _F(getMaterialProperty<RealTensorValue>("deformationGradient")),
    _invFtr(getMaterialPropertyOld<RealTensorValue>("invFtr")),
    _invFtr_old(getMaterialPropertyOld<RealTensorValue>("invFtr")),
    _kappa(getParam<Real>("kappa")),
    _alpha_f(getParam<Real>("alpha"))
{

   for (unsigned int i = 0; i < _numComp; ++i) 
        {
         _disp_var[i] = coupled("displacements",i);
        }
  
    if(_var.number() == (_disp_var[0])) component=0;
    if(_var.number() == (_disp_var[1]) && _numComp >= 2) component=1;
    if(_var.number() == (_disp_var[2]) && _numComp == 3)  component=2;

   for (int i = 0; i < _numComp; ++i) zeros(i) = 0.0;
   
   _J_lin = new Real * [64];
   
   for (int j=0; j< 64; ++j)_J_lin[j]=new Real[64];
   
   _Pvol_Lin = new RealTensorValue * [64];
   
   for (int j=0; j< 64; ++j) _Pvol_Lin[j]=new RealTensorValue[64];
    
}

Real HHTNearlyIncompressibility::computeQpResidual() {
    
  //  RealTensorValue _Cinv  = _F[_qp].inverse() * _invFtr[_qp];
    RealTensorValue  P_vol =  _kappa * (1.0 - _alpha_f) * (_J[_qp]-1.0) * _J[_qp] * _invFtr[_qp]
                             + _alpha_f * _kappa * (_J_old[_qp]-1.0) * _J_old[_qp] * _invFtr_old[_qp]
    ;
    return   ( P_vol(component , 0) * _grad_test[_i][_qp] (0)+
               P_vol(component , 1) * _grad_test[_i][_qp] (1)+
               P_vol(component , 2) * _grad_test[_i][_qp] (2));
}

Real HHTNearlyIncompressibility::computeQpJacobian()
{     

    if (_i==0 && _j==0 && _qp==0)
      {
          H = RealTensorValue(zeros, zeros, zeros);
          
          for (int j=0; j< _phi.size() ; ++j)
              for (int qp=0; qp< _qrule->n_points() ; ++qp)
              {
                  for (unsigned k = 0; k < _numComp; ++k)
                  {
                      H(component, k) = _grad_phi[j][qp](k);
                  }
                  _Pvol_Lin[j][qp] = _kappa *  (1.0 - _alpha_f) * ( 2.0 * _J[qp] * _J[qp] *_invFtr[qp].contract(H) * _invFtr[qp]
                                               -1.0 * _J[qp]  * _invFtr[qp].contract(H) * _invFtr[qp] +
                                               -1.0 * _J[qp]  * _J[qp] * _invFtr[qp] * H.transpose() * _invFtr[qp] +
                                                1.0 * _J[qp]  * _invFtr[qp] * H.transpose() * _invFtr[qp] );
              }
        }
        
     
          
         return  _Pvol_Lin[_j][_qp](component,0) * _grad_test[_i][_qp](0)+
                 _Pvol_Lin[_j][_qp](component,1) * _grad_test[_i][_qp](1)+
                 _Pvol_Lin[_j][_qp](component,2) * _grad_test[_i][_qp](2);
    


}  
      
       



Real HHTNearlyIncompressibility::computeQpOffDiagJacobian(unsigned int jvar) {
   
    if (jvar==_disp_var[0] || jvar==_disp_var[1] || jvar ==_disp_var[2])
    {
        
        unsigned int disp_index = 99;

        if (jvar == _disp_var[0])
            disp_index = 0;
        
        else if (jvar == _disp_var[1])
            disp_index = 1;
        
        else if (jvar == _disp_var[2])
            disp_index = 2;
    
        if (_i == 0 && _j == 0 && _qp == 0)
        {
          H = RealTensorValue(zeros, zeros, zeros);

          // for all the components of the trial functions:
          for (int j = 0; j < _phi.size(); ++j)

            // for all the  quadrature points:
            for (int qp = 0; qp < _qrule->n_points(); ++qp)
            {

              // in H we store the gradient of trial function: in this way j becomes
              // _jvar for H
              for (unsigned k = 0; k < _numComp; ++k)
              {
                H(disp_index, k) = _grad_phi[j][qp](k);
              }

                        _Pvol_Lin[j][qp] = _kappa * (1.0 - _alpha_f) * ( 2.0 * _J[qp] * _J[qp] *_invFtr[qp].contract(H) * _invFtr[qp]
                                                     -1.0 * _J[qp] * _J[qp] * _invFtr[qp] * H.transpose() * _invFtr[qp] +
                                                     -1.0 * _J[qp] * _invFtr[qp].contract(H) * _invFtr[qp] +
                                                      1.0 * _J[qp] * _invFtr[qp] * H.transpose() * _invFtr[qp] );
             }
        }
            
         
              
             return  _Pvol_Lin[_j][_qp](component,0) * _grad_test[_i][_qp](0)+
                     _Pvol_Lin[_j][_qp](component,1) * _grad_test[_i][_qp](1)+
                     _Pvol_Lin[_j][_qp](component,2) * _grad_test[_i][_qp](2);
    }
    
    else
        return 0.0;

 
}

