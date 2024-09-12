[Mesh]
[file]
  type = FileMeshGenerator
  file = hydrocoin.e
  block_id = '1 6 4 7 2'
  boundary_id = '6 7'
  boundary_name = 'inflow outflow'
[]

[inclusionsList]
 input = file
 type = InclusionsMeshModifier
 fn = 2
 fx_string = '0.0 0.0'
 fy_string = '-6500.0 500.0'
 fa1_string = '79.695153531233970 -45.0'
 fd1_string = '20000000.0 200000000.0'
 fd2_string = '14.758048690536125 7.071067811865476'
[]

[inclusionRefinement]
 input = inclusionsList
 type  = InclusionRefinement
 inclusions_list = inclusionsList
 doBoundaryRefinement = true
 refinements = '6 0'  # amr uniforme
 outputFileName = test.e
[]
[]

[Variables]
[./pressure] order=FIRST  family=LAGRANGE [../]
[]

[Kernels]
[./diffusion] 
type = PermeabilityDiffusion 
variable = pressure 
[../]
[]


[Materials]
[./flowAndTransport] 
type = FlowAndTransport 
k = 1 kFrac = 1e2
inclusions_list = inclusionsList 
[../]
[]

# observe that with the second BCs the stiffness matrix is SDP and we can use choleski factorization
[BCs]
[./inflowBC]  type = FunctionDirichletBC variable = pressure function = y  boundary = 1  [../]
[]
 

[Preconditioning]
[./prec] type = SMP full = true ksp_norm = default [../]
[]

[Executioner]

 type=Steady
 solve_type= LINEAR
 line_search = none
#petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
#petsc_options_value='  preonly   lu       NONZERO               mumps '
 
 petsc_options_iname = '-ksp_monitor_singular_value'
 petsc_options_value = '::ascii_matlab'

[./Quadrature]
 type=GRID
 order = FIFTH
[../]

[]


[Outputs]
 file_base  = DiffusionOut_4
 exodus     = true
 perf_graph = true
[]

