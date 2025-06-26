#ifndef FIBERSGEOMETRY_H
#define FIBERSGEOMETRY_H

#include "FibersMaterial.h"

class FibersGeometry;

template <>
InputParameters validParams<FibersGeometry>();

/**
 * Material for providing an interface to fibre direction
*
 * Material properties are:
 *  - the local coordinate system's basis vectors \f$\hat{e}_f\f$,
 * \f$\hat{e}_n\f$, \f$\hat{e}_s\f$
 *    (see notation in \ref Holzapfel2009 "Holzapfel, 2009, Figure 1")
 *  - an appropriate rotation matrix \f$\mathbf{R}=(\hat{e}_f, \hat{e}_n,
 * \hat{e}_s)\f$, i.e.
 *    containing the unit vectors column-wise
 * For debugging purposes, the rotation matrix (and respective coordinate
 * system) can also be given externally in the input file.
 *
 * The local fibre coordinate system is constructed in an analogous
 * fashion to the description in \ref Potse2006
 * with the difference that we do not average over neighbouring elements
 * for getting a smoothed thickness parameter e.
 * Instead, we directly use the result from a CardiacThicknessParameterAux
 * kernel.
 * For further details consult the publication and see into the
 * CardiacThicknessParameterAux and VolumeNearestNodeAux kernels.
 */
class FibersGeometry : public FibersMaterial
{
public:
  FibersGeometry(const InputParameters &parameters);

protected:
  virtual void computeQpProperties();

private:
  /// Computes \f$\frac{\vec{v}}{|\vec{v}|}\f$.
  inline RealVectorValue
  VectorNormalize(const RealVectorValue &v)
  {
    return v / v.norm();
  }

  const VariableValue &_e;
  const VariableGradient &_grad_e;

  RealVectorValue E_fiber;
  RealVectorValue E_sheet;
  RealVectorValue E_normal;
};

#endif

