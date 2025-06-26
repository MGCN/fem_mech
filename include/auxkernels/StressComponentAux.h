//* This file is part of the MOOSE framework
//
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

// Forward Declarations

/**
 * Computes h_min / |u|
 */
class StressComponentAux : public AuxKernel
{
public:
  static InputParameters validParams();

  StressComponentAux(const InputParameters & parameters);

  virtual ~StressComponentAux() {}

protected:
  virtual Real computeValue();

  // Velocity gradients
  //FEProblem * _fe_problem_ptr;
  const VariableGradient & _grad_velocity_x;
  const VariableGradient & _grad_velocity_y;
  const VariableGradient & _grad_velocity_z;
  double _mu;
  const bool _use_normal;

  /// Will hold 0, 1, or 2 corresponding to x, y, or z.
  const int _component;
  /// normals at quadrature points
  const MooseArray<Point> & _normals;
  //int _boundary_tags;
};
