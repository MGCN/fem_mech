

#ifndef StressAux_H
#define StressAux_H

#include "AuxKernel.h"
#include "ElasticMaterial.h"
class StressAux;

template <>
InputParameters validParams<StressAux>();

class StressAux : public AuxKernel

{
public:
  StressAux(const InputParameters &parameters);

  ~StressAux(){};

public:
    

protected:

   virtual Real computeValue();

  FEProblem * _fe_problem_ptr;
  /// normals at quadrature points
   const MooseArray<Point> & _normals;
   
   const MaterialProperty<RealTensorValue> &_stress;
 
   int _boundary_tags;
};
#endif

