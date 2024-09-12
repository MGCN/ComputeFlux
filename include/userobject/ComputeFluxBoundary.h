//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "GeneralUserObject.h"
#include "Function.h"
#include "SystemBase.h"

// LibMesh includes
#include "libmesh/linear_implicit_system.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/mesh.h"
#include "libmesh/distributed_mesh.h"

// Parrot includes
// #include "FlowAndTransport.h"
// #include "ParrotSolver.h"


class ComputeFluxBoundary;

template <>
InputParameters validParams<ComputeFluxBoundary>();

class ComputeFluxBoundary : public GeneralUserObject
{
public:
  static InputParameters validParams();

  ComputeFluxBoundary(const InputParameters & params);

  virtual void initialize ()                   override;
  virtual void execute    ()                   override {};
  virtual void finalize   ()                   override {};
  virtual void threadJoin (const UserObject &) override {};

  protected:

  void init();
  void assembleMatrix();
  void assembleFracture();
  void setDirichletMatrix();
  void setDirichletFracture();
  void solve();
  // void write();
  void del();
  Real ComputeMaterialProprties(const Elem *elem);
  Real ComputeMaterialProprtiesFracture(const Elem *elem);

  bool _material;
  std::string  _material_name;
  std::string  _var_name_M;
  std::string  _sys_name_M;
  std::string  _mat_domain_name_M0;
  std::string  _mat_boundary_name_M0;
  std::string  _vec_flux_name_M0;
  std::string  _vec_rhs_name_M0;
  std::string  _vec_boundary_name_M0;
  std::string  _vec_dir_M;


  std::string  _var_name_F;
  std::string  _sys_name_F;

  std::string  _mat_domain_name_F0;
  std::string  _mat_boundary_name_F0;
  std::string  _vec_flux_name_F0;
  std::string  _vec_rhs_name_F0;
  std::string  _vec_boundary_name_F0;
  std::string  _vec_dir_F;

  QBase const * const & _qrule;

  int _block_0_id;
  int _block_1_id;
  int _block_2_id;

  EquationSystems       & _equationSystemsM;
  MeshBase         * _meshBase;
  MeshBase         * _meshBaseF;
  LinearImplicitSystem  * _linearImplicitSystemM;

  LinearImplicitSystem  * _linearImplicitSystemF;


  unsigned int _v_varM;

  unsigned int _v_varF;

  MaterialBase     * _materialBase;
  // FlowAndTransport * _flowAndTransport;

  unsigned int _dim_M;

  unsigned int _dim_F;

  DofMap * _dof_map;

  DofMap * _dof_mapF;

  std::vector<BoundaryName>     _dirichletNodesetNames_M;
  std::vector<BoundaryID>       _dirichletSidesetIds_M;
  std::vector<Function const *> _dirichletFunctions_M;
  std::vector<FunctionName>     _dirichletFunctionNames_M;


  std::vector<BoundaryName>     _dirichletNodesetNames_F;
  std::vector<BoundaryID>       _dirichletSidesetIds_F;
  std::vector<Function const *> _dirichletFunctions_F;
  std::vector<FunctionName>     _dirichletFunctionNames_F;


  std::vector<BoundaryName>     _neumannSidesetNames_M;
  std::vector<FunctionName>     _neumannFunctionNames_M;

  std::vector<BoundaryName>     _neumannSidesetNames_F;
  std::vector<FunctionName>     _neumannFunctionNames_F;

  std::vector<BoundaryID>       _neumannSidesetIds_M;
  std::vector<BoundaryID>       _neumannSidesetIds_F;

  std::vector<Function const *> _neumannFunctions_M;
  std::vector<Function const *> _neumannFunctions_F;


  bool const _has_exodus_file;
  bool const _has_variable_name_1;
  bool const _has_variable_name_2;
  bool const _has_block_2_id;

  std::string _exodus_filename;
  std::string _variable_name_1;
  std::string _variable_name_2;
  std::vector<dof_id_type>  _dirIds;
  std::vector<int> _vector_p;
  std::vector<Real> _vector_value_m;
  std::vector<Real> _vector_value_f;
  std::string _multiapp_name;
  

};

