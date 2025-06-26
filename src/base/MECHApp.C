//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MECHApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

// So we can register objects from the fluid_properties module.


template <>
InputParameters
validParams<MECHApp>()
{
    InputParameters params = validParams<MooseApp>();
    
    return params;
}

registerKnownLabel("MECHApp");

MECHApp::MECHApp(const InputParameters & parameters) : MooseApp(parameters)
{
    Moose::registerObjects(_factory);
    MECHApp::registerObjects(_factory);
    
    Moose::associateSyntax(_syntax, _action_factory);
    MECHApp::associateSyntax(_syntax, _action_factory);
    
}

MECHApp::~MECHApp() {}

extern "C" void
MECHApp_registerApps()
{
    MECHApp::registerApps();
}
void
MECHApp::registerApps()
{
    registerApp(MECHApp);
}

// External entry point for dynamic object registration
extern "C" void
MECHApp__registerObjects(Factory & factory)
{
    MECHApp::registerObjects(factory);
}
void
MECHApp::registerObjects(Factory & factory)
{
    Registry::registerObjectsTo(factory, {"MECHApp"});
}

// External entry point for dynamic syntax association
extern "C" void
MECHApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
    MECHApp::associateSyntax(syntax, action_factory);
}
void
MECHApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
      Registry::registerActionsTo(action_factory, {"MECHApp"});
}

