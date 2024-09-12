#include "parrot2App.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
parrot2App::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy DirichletBC, that is, set DirichletBC default for preset = true
  //params.set<bool>("use_legacy_dirichlet_bc") = false;

  return params;
}

parrot2App::parrot2App(InputParameters parameters) : MooseApp(parameters)
{
  parrot2App::registerAll(_factory, _action_factory, _syntax);
}

parrot2App::~parrot2App() {}

void
parrot2App::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"parrot2App"});
  Registry::registerActionsTo(af, {"parrot2App"});

  /* register custom execute flags, action syntax, etc. here */
}

void
parrot2App::registerApps()
{
  registerApp(parrot2App);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
parrot2App__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  parrot2App::registerAll(f, af, s);
}
extern "C" void
parrot2App__registerApps()
{
  parrot2App::registerApps();
}
