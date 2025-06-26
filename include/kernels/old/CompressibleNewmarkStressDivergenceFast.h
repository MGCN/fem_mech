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
/****************************************************************/



#ifndef CompressibleNewmarkStressDivergenceFast_H
#define CompressibleNewmarkStressDivergenceFast_H

#include "Kernel.h"
#include "ElasticMaterial.h"

class CompressibleNewmarkStressDivergenceFast;

template <>
InputParameters validParams<CompressibleNewmarkStressDivergenceFast>();

class CompressibleNewmarkStressDivergenceFast : public Kernel
{
public:
  CompressibleNewmarkStressDivergenceFast(
                                InputParameters const & parameters);



protected:

    virtual void computeResidual();
    virtual void computeResidual2D();
    virtual void computeResidual3D();
    
    virtual  void computeJacobian(){};
    virtual  void computeJacobian2D();
    virtual  void computeJacobian3D();
    
    virtual  void computeOffDiagJacobian(unsigned int );
    
    virtual Real computeQpResidual(){return 0.0;};

    virtual Real computeQpJacobian(){return 0.0;};

    virtual Real computeQpOffDiagJacobian(unsigned int){return 0.0;};
    
    
    unsigned const _dim;

    unsigned const _disp_x_var;
    unsigned const _disp_y_var;
    unsigned const _disp_z_var;

    
    RealVectorValue zeros;

    RealTensorValue H;

    
    const MaterialProperty<RealTensorValue> &_P;
    const MaterialProperty<RealTensorValue> &_P_old;
    const MaterialProperty<RealTensorValue> &_P_older;
    
    RealTensorValue ***_stress_lin;

    
    const MaterialProperty<ElasticMaterial *> &_materialPointer;



  
};

#endif // CompressibleNewmarkStressDivergenceFast_H
