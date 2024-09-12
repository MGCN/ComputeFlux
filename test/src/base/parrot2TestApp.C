//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "parrot2TestApp.h"
#include "parrot2App.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
parrot2TestApp::validParams()
{
  InputParameters params = parrot2App::validParams();
  return params;
}

parrot2TestApp::parrot2TestApp(InputParameters parameters) : MooseApp(parameters)
{
  parrot2TestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

parrot2TestApp::~parrot2TestApp() {}

void
parrot2TestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  parrot2App::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"parrot2TestApp"});
    Registry::registerActionsTo(af, {"parrot2TestApp"});
  }
}

void
parrot2TestApp::registerApps()
{
  registerApp(parrot2App);
  registerApp(parrot2TestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
parrot2TestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  parrot2TestApp::registerAll(f, af, s);
}
extern "C" void
parrot2TestApp__registerApps()
{
  parrot2TestApp::registerApps();
}
