#ifndef MECHAPP_H
#define MECHAPP_H

#include "MooseApp.h"

class MECHApp;

template <>
InputParameters validParams<MECHApp>();

class MECHApp : public MooseApp
{
public:
  MECHApp(const InputParameters &parameters);
  virtual ~MECHApp();

  static void registerApps();
  static void registerObjects(Factory &factory);
  static void associateSyntax(Syntax &syntax, ActionFactory &action_factory);
};

#endif /* MECHAPP_H */
