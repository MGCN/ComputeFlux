[Problem]
# solve = false
# skip_nl_system_check = true
[]

# [Problem]
# type = ParrotProblem3
# use_AFC = true
# operator_userobject = storeOperatorsUO
# solver_type = 1
# []


[Mesh]
[gmg]
 type = GeneratedMeshGenerator 
 dim = 2
 nx = 80
 ny = 80
 #type = FileMeshGenerator
 #file = adapted_640_4_0000_mesh.xdr
[ ]

[./subdomains_0]
input = gmg
type = SubdomainBoundingBoxGenerator
bottom_left = '0.0 0.00001 0.0'
top_right =   '0.650000 1.0001 0.0'
block_id = 0
# block_name = up_block_0
[../]


[./subdomains_1]
input = subdomains_0
type = SubdomainBoundingBoxGenerator
bottom_left = '0.6500001 0.00001 0.0'
top_right =   '1.000 1.0001 0.0'
block_id = 1
#block_name = up_block_1
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
[]

[inclusionRefinement]
input = inclusionsList
type  = InclusionRefinement
inclusions_list = inclusionsList
doBoundaryRefinement = false
refinements = '8 0'
#outputFileName = test.e
[]
[]


[Variables]
[pressure] 
[../]
[]

[AuxVariables]
[flux_1]
[../]
[flux_2]
[../]
[u_aux]
[../]

[]

# [Kernels]
# [upwind] type = Advection variable = CM [../]
# [./time] type = PorosityTimeDerivative variable = CM lumping = true [../]
# []

# [Materials]
# [./flowAndTransport] type = FlowAndTransport k = 1.0 kFrac = 1.0e4 phi = 1.0 phiFrac = 1.0 
# pressure=pressure2 inclusions_list = inclusionsList [../]
# []


# [BCs]
# [./u_injection_left] type = DirichletBC boundary = left  variable = CM value=1 [../]
# []

[Kernels]
[permeabilityDiffusion] type = PermeabilityDiffusion variable = pressure [../]
[]

[Materials]
[flowAndTransport]
type = FlowAndTransport 
#block = 0
k = 1.0
phi = 0.2 
kFrac = '1.0e4'
phiFrac = 0.4
pressure = pressure
inclusions_list = inclusionsList # inclusionsListRefinement
[]
[]

[BCs]
[dirinflow] type = DirichletBC variable = pressure value = 1 boundary = 1 []
[diroutflow] type = NeumannBC variable = pressure value = 1 boundary = 3 []
[]

[Executioner]

 type = Steady
 solve_type= LINEAR
 line_search = none
 petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 petsc_options_value='  preonly   lu       NONZERO               mumps '
 petsc_options = '-ksp_monitor_singular_value '

 [Quadrature] type = GRID   order = ELEVENTH []

[]

[Outputs]
file_base = AdvectionOut_flow
exodus = true
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
# nemesis_file = outputUO.e
exodus_file  = outputUO2.e
variable_name = 'pressure2'
[../]

[./assF]
type = ComputeFlux
execute_on = 'timestep_end'
block_0_id='0'
block_1_id='1'
material_name = flowAndTransport
dirichlet_nodeset_names  = ' 1 '
dirichlet_function_names = ' gamma1'
neumann_sideset_names  = ' 3'
neumann_function_names = 'bottom2 '
variable_name_1='flux_1'
variable_name_2='flux_2'
#exodus_file=flux.e
#mesh=porous_matrix_out.e
[../]

[./storeOperatorsUO]
type = StoreOperators
[../]

[]



[Functions]
 [gamma1]  type = ParsedFunction value = 1.0 []
 [bottom2] type = ParsedFunction value = 1.0 []
[]



