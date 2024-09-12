[Problem]
[]



[Mesh]
[gmg]
 type = FileMeshGenerator
 file = back_80.xda
[../]



[./subdomains_0]
input = gmg
type = SubdomainBoundingBoxGenerator
bottom_left = '0.0   0.0  0.0'
top_right =   '1.0   0.70  0.0'
block_id = 0
[../]


[./subdomains_1]
input = subdomains_0
type = SubdomainBoundingBoxGenerator
bottom_left = '0.0   0.70 0.0'
top_right =   '1.0001 1.0001 0.0'
block_id = 1
[../]



[inclusionsList]
input = subdomains_1
type = InclusionsMeshModifier
fn = 6
fx_string = '0.5 0.5 0.75 0.75 0.625 0.625'
fy_string = '0.5 0.5 0.75 0.75 0.625 0.625'
fd1_string = '1.0 1.0 0.5 0.5 0.25 0.25'
fd2_string = '1.0e-4 1.0e-4 1.0e-4 1.0e-4 1.0e-4 1.0e-4'
fa1_string = '0.0 90.0 0.0 90.0 0.0 90.0'
[../]

[inclusionRefinement]
input = inclusionsList
type  = InclusionRefinement
inclusions_list = inclusionsList
doBoundaryRefinement = false
refinements = '0 0'
[../]
[]


[Variables]
[pressure] 
[../]
[]

[AuxVariables]
[u_aux]
[../]
[]


[Kernels]
[permeabilityDiffusion] 
type = PermeabilityDiffusion 
variable = pressure 
[../]
[]

[Materials]
[flowAndTransport]
type = FlowAndTransport 
k = 1.0
phi = 0.2 
kFrac = '1.0e4'
phiFrac = 0.4
pressure = pressure
inclusions_list = inclusionsList # inclusionsListRefinement
[../]
[]

[BCs]
[dirinflow] 
type = DirichletBC 
variable = pressure 
value = 1 
boundary = 1 
[../]
[diroutflow] 
type = NeumannBC 
variable = pressure 
value = 1 
boundary = 3 
[../]
[]

[Executioner]
type = Steady
solve_type= LINEAR
line_search = none
petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
petsc_options_value='  preonly   lu       NONZERO               mumps '
petsc_options = '-ksp_monitor_singular_value '

[Quadrature] 
type = GRID   
order = ELEVENTH 
[../]
[]

[Outputs]
file_base = AdvectionOut_flow_4000
nemesis = true
[]


[UserObjects]
[./my]
type = SolveDiffusion
execute_on = initial
material_name = flowAndTransport
dirichlet_nodeset_names  = '1 '
dirichlet_function_names = ' gamma1 '
neumann_sideset_names  = ' 3 '
neumann_function_names = 'bottom2 '
solver_type = 3
variable_name = 'pressure2'
[../]

[./assF]
type = ComputeFluxBC
execute_on = 'timestep_end'
material_name = flowAndTransport
dirichlet_nodeset_names  = ' 1 '
dirichlet_function_names = ' gamma1'
neumann_sideset_names  = ' 3'
neumann_function_names = 'bottom2 '
[../]

[./storeOperatorsUO]
type = StoreOperators
[../]

[]



[Functions]
 [gamma1]  type = ParsedFunction value = 1.0 [../]
 [bottom2] type = ParsedFunction value = 1.0 [../]
[]



