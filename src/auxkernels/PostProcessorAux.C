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
/* Auxiliary Kernel to visualize the fibers.                    */
/****************************************************************/


/* just for P1 */

#include "PostProcessorAux.h"
registerMooseObject("MECHApp", PostProcessorAux);
template <>
InputParameters
validParams<PostProcessorAux>()
{
  InputParameters params = validParams<AuxKernel>();

  MooseEnum component(
      "E_xx E_yy E_zz E_yz E_xz E_xy S_xx S_yy S_zz S_yz S_xz S_xy J VonMises");

  params.addRequiredParam<MooseEnum>("component", component,
                                     "The desired value.");
  return params;
}

PostProcessorAux::PostProcessorAux(const InputParameters &parameters)
    : AuxKernel(parameters), _component(getParam<MooseEnum>("component")),
      _stress(getMaterialProperty<RealTensorValue>("CauchyStress")),
      _F(getMaterialProperty<RealTensorValue>("deformationGradient")),
     _invFtr(getMaterialProperty<RealTensorValue>("invFtr")),
     _kappa(getParam<Real>("kappa"))
{
}

Real
PostProcessorAux::computeValue()
{
  RealTensorValue F = _F[_qp];

  RealTensorValue Ftr = F.transpose();
  
  Real J = F.det();

  RealTensorValue P_vol = _kappa * (J - 1.0) * J * _invFtr[_qp];

  RealTensorValue sigma_vol= P_vol * _F[_qp].transpose()/J;

  RealTensorValue C = F.transpose() * F;

  RealTensorValue E = C;

  for (int i = 0; i < 3; ++i)
  {
    E(i, i) -= 1.0;
  }

  E = 0.5 * E;

  RealTensorValue S = _stress[_qp];

  RealTensorValue TEST = S - S.transpose();
  Real sum = 0.0;
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      sum += std::fabs(TEST(i, j));
    }
  }

  switch (_component)
  {
    case 0:
    {
      return E(0, 0);
    }
    case 1:
    {
      return E(1, 1);
    }
    case 2:
    {
      return E(2, 2);
    }
    case 3:
    {
      return E(1, 2);
    }
    case 4:
    {
      return E(0, 2);
    }
    case 5:
    {
      return E(0, 1);
    }
    case 6:
    {
      return S(0, 0);
    }
    case 7:
    {
      return S(1, 1);
    }
    case 8:
    {
      return S(2, 2);
    }
    case 9:
    {
      return S(1, 2);
    }
    case 10:
    {
      return S(0, 2);
    }
    case 11:
    {
      return S(0, 1);
    }
    case 12:
    {
      return J;
    }
    case 13:
    {
      Real VM = std::sqrt(
          0.5 *
          (std::pow(S(0, 0) - S(1, 1), 2.0) + std::pow(S(2, 2) - S(1, 1), 2.0) +
           std::pow(S(0, 0) - S(2, 2), 2.0) +
           6.0 * (S(0, 1) * S(0, 1) + S(2, 1) * S(2, 1) + S(0, 2) * S(0, 2))));
 
//      std::cout<<"ciao"<<'\n';
      return VM;
    }
    default:
    {
      mooseWarning("Unexpected component in PostProcessorAux");
      return 0.0;
    }
  }
}
