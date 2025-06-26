
#ifndef GrowthTrace_H
#define GrowthTrace_H

#include "AuxKernel.h"
#include "ElasticMaterial.h"
#include "GrowthMaterial.h"

class GrowthTrace;

template <>
InputParameters validParams<GrowthTrace>();

class GrowthTrace : public AuxKernel

{
public:
  GrowthTrace(const InputParameters &parameters);

  ~GrowthTrace(){};

public:
    

protected:
 
   virtual Real computeValue();
   
   const MaterialProperty<RealTensorValue> &_Se;

   const MaterialProperty<RealTensorValue> &_Ce;
 
   // int _boundary_tags;
   
};
#endif
