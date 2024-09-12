[Problem]
# solve = false
[]

[Mesh]
[gmg]
 type = FileMeshGenerator
 # file = mfd_complex_fine.msh
 file = mfd_complex.msh
[]
second_order = true
[nodeset1]
 input = gmg 
 type = BoundingBoxNodeSetGenerator
 new_boundary = top
 top_right =   ' 1.1 1.1 0.0'
 bottom_left = '-0.1 0.9 0.0'
[]

[nodeset2]
 input = nodeset1 
 type = BoundingBoxNodeSetGenerator
 new_boundary = bottom
 top_right =   ' 1.1  0.1 0.0'
 bottom_left = '-0.1 -0.1 0.0'
[]


[inclusionsList]
 input = nodeset2
 type = InclusionsMeshModifier
fn = 10
fx_string = ' 0.135 0.15 0.3 0.275 0.74986 0.77486 0.725 0.575 0.85 0.275'
fy_string = ' 0.2392 0.205 0.36 0.70835 0.50046 0.20131 0.32375 0.84285 0.88645 0.9045'
fd1_string = ' 0.39234 0.24413 0.61774 0.48594 0.69499 0.16418 0.27415 0.51827 0.24523 0.28479'
fd2_string = ' 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001'
fa1_string = ' 115.6768 145.008 119.0546 120.9617 106.7009 155.7724 155.7723 150.2592 144.6443 208.6169'
[../]

[]

[Variables]
[pressure] order = second []
# [pressure] order = THIRD family = HIERARCHIC []
[]

[Kernels]
[./diffusion] type = PermeabilityDiffusion variable = pressure [../]
[]

[Materials]
[./flowAndTransport]
type = FlowAndTransport k = 1
kFrac = '1e4 1e4 1e4 1e-4 1e-4 1e4 1e4 1e4 1e4 1e4'
inclusions_list = inclusionsList
operation_type = harmonic
[../]
[]

[BCs]
[./bot] type = DirichletBC variable = pressure value = 4.0 boundary = top    [../]
[./top] type = DirichletBC variable = pressure value = 1.0 boundary = bottom [../]
[]

[Executioner]
 type = Steady
 solve_type = LINEAR
 # petsc_options_iname = '-pc_type -pc_hypre_type '
 # petsc_options_value = ' hypre    boomeramg     '
 petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 petsc_options_value='  preonly   lu       NONZERO               mumps '


[Quadrature] type = GRID order = SEVENTH []

[]

[Outputs]
 file_base = DiffusionOut_R
 nemesis = true
[]
