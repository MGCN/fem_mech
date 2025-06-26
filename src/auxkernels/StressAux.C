#include "StressAux.h"
#include "FEProblem.h"
#include "Assembly.h"
#include "libmesh/quadrature_gauss.h"
#include <limits>
double Inf = std::numeric_limits<double>::infinity();
registerMooseObject("MECHApp", StressAux);

template <>


InputParameters
validParams<StressAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addParam<int>("boundary_ids", "The list of boudary ID");
  return params;
}

StressAux::StressAux(const InputParameters &parameters)
    :  AuxKernel(parameters),
      _fe_problem_ptr(getParam<FEProblem *>("_fe_problem")),
      _normals(_assembly.normals()),
      _stress(getMaterialProperty<RealTensorValue>("firstPiolaKirchhofTensor"))
      
{
      if (isParamValid("boundary_ids")){
        _boundary_tags=getParam<int>("boundary_ids");
       }
}


Real
StressAux::computeValue()
{
  Real result=0;  
    
    for(uint sideno = 0; sideno < _current_elem->n_sides(); ++sideno)
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
               
               fe_elem_face->reinit(_current_elem, sideno);
               
               if(qface_normals[_qp](0) || qface_normals[_qp](1) || qface_normals[_qp](2))
                  {
                     RealVectorValue temp = _stress[_qp] * qface_normals[_qp];
                     Real result_2 = qface_normals[_qp] * temp;
                     if (std::abs(result_2)<1.e10 && std::abs(result_2)!=Inf && std::abs(result_2)!=NAN)
                     result = qface_normals[_qp] * temp;

                 }
             }
         }
     }
    
  
 // RealVectorValue temp = _stress[_qp] * _normals[_qp];

 // result = _normals[_qp] * temp;
  
  return result;

}



