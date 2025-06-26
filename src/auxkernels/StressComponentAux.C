//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "StressComponentAux.h"
#include "MooseMesh.h"
#include "Assembly.h"
registerMooseObject("MECHApp", StressComponentAux);

InputParameters
StressComponentAux::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addClassDescription("This class computes the wall shear stressbased on "
                             "pressure and velocity for incompressible Navier-Stokes");
  params.addCoupledVar("velocity_x", "The velocity component");
  params.addCoupledVar("velocity_y", "The velocity component");
  params.addCoupledVar("velocity_z", "The velocity component");
  //params.addParam<int>("boundary_ids", "The list of boudary ID");
  // params.addRangeCheckedParam<unsigned int>("comp", 0, "0<=comp<=2", "The component");
  MooseEnum component("x y z normal");
  params.addRequiredParam<MooseEnum>("component", component, "The desired component of flux.");
  params.addParam<double>("mu", "The viscosity");

  return params;
}

StressComponentAux::StressComponentAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    //_use_normal(getParam<MooseEnum>("component") == "normal"),
    //_fe_problem_ptr(getParam<FEProblem *>("_fe_problem")),
    _grad_velocity_x(isCoupled("velocity_x") ? coupledGradient("velocity_x") : _grad_zero),
    _grad_velocity_y(isCoupled("velocity_y") ? coupledGradient("velocity_y") : _grad_zero),
    _grad_velocity_z(isCoupled("velocity_z") ? coupledGradient("velocity_z") : _grad_zero),
    _mu(getParam<double>("mu")),
    _use_normal(getParam<MooseEnum>("component") == "normal"),
    _component(getParam<MooseEnum>("component")),
    _normals(_assembly.normals())
{
	
      if (_use_normal && !isParamValid("boundary")){
        paramError("boundary", "A boundary must be provided if using the normal component!");
       }
}

Real
StressComponentAux::computeValue()
{
	
  //std::cout<<"Compute WSS"<<std::endl;

  RealTensorValue v(_grad_velocity_x[_qp], _grad_velocity_y[_qp], _grad_velocity_z[_qp]);

  //std::cout<<"velocity"<<v<<std::endl;

  Real value_wss = 0;  

  RealTensorValue sigma = 0.5 * _mu * (v + v.transpose());
    
/*  for(uint sideno = 0; sideno < _current_elem->n_sides(); ++sideno)
    {
        int t = _boundary_tags;
        {
           if (_mesh.getMesh().get_boundary_info().has_boundary_id(_current_elem,sideno,t))
           { 
              

              FEType fe_type = _fe_problem_ptr->getStandardVariable(0,"disp_x").feType();
              
              int dim = _mesh.getMesh().mesh_dimension();
               
               std::unique_ptr<FEBase> fe_elem_face(FEBase::build(dim, fe_type));

               
               QGauss qface(dim-1, fe_type.default_quadrature_order());
               
               fe_elem_face->attach_quadrature_rule(&qface);

               const std::vector<Point> & qface_normals = fe_elem_face->get_normals();
               
               //fe_elem_face->reinit(_current_elem, sideno);

	       std::cout<<"Compute WSS 2"<<qface_normals[_qp]<<std::endl;
       
       	       if(qface_normals[_qp](0) || qface_normals[_qp](1) || qface_normals[_qp](2))
                 {
                     RealVectorValue traction = sigma * qface_normals[_qp];
                     
		     Real traction_n = qface_normals[_qp] * traction;

		     RealVectorValue WSS = traction - traction_n * qface_normals[_qp];
  
  		     //if (std::abs(result_2)<1.e10 && std::abs(result_2)!=Inf && std::abs(result_2)!=NAN)
  
	             //result = qface_normals[_qp] * temp;

		     value_wss =std::sqrt((std::pow(WSS(0),2)) + (std::pow(WSS(1),2)) + (std::pow(WSS(2),2)));
		     //
		     std::cout<<"value_wss"<<value_wss<<std::endl;

                 }
              }
           }
        }
    
  */  
  

   //std::cout<<"normals"<<_normals[_qp]<<std::endl;

   RealVectorValue sigma_n = sigma * _normals[_qp];

   //std::cout<<"sigma_n"<<sigma_n(0)<<std::endl;

   RealVectorValue WSS  = sigma_n - (sigma_n * _normals[_qp] ) * _normals[_qp];

   //std::cout<<"WSS"<<WSS(0)<<std::endl;

   value_wss = std::sqrt((std::pow(WSS(0),2)) + (std::pow(WSS(1),2)) + (std::pow(WSS(2),2)));

  
   //std::cout<<"WSS"<< value_wss<<std::endl;
  
   return value_wss;
}
