#ifndef PostProcessorAux_H
#define PostProcessorAux_H

#include "AuxKernel.h"

class PostProcessorAux;

template <>
InputParameters validParams<PostProcessorAux>();

class PostProcessorAux : public AuxKernel
{
public:
  PostProcessorAux(const InputParameters &parameters);

  ~PostProcessorAux(){};

protected:
  virtual Real computeValue();

  unsigned _component;

  const MaterialProperty<RealTensorValue> &_stress;
  const MaterialProperty<RealTensorValue> &_F;
  const MaterialProperty<RealTensorValue> &_invFtr;
  Real _kappa;
};
#endif
