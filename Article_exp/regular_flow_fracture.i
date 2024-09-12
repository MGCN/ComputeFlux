[Mesh]
[gmg]
type = FileMeshGenerator
file = fracture_network_out.e
[]


[./subdomains_0]
input = gmg
type = SubdomainBoundingBoxGenerator
bottom_left = '0.0   0.0  0.0'
top_right =   '0.7   1.0 0.0'
block_id = 1
[../]


[./subdomains_1]
input = subdomains_0
type = SubdomainBoundingBoxGenerator
bottom_left = '0.70001 0.0 0.0'
top_right =   '1.0001 1.0001 0.0'
block_id = 2
[../]


[]

[Problem]
type = FEProblem
solve=false
[]


[Variables]
[pressure_F] 
 [../]
[]

[AuxVariables]
[flux_1]
[../]
[flux_2]
[../]
[aux_sol_F]
[../]

[]
 

[Kernels]
[permeabilityDiffusion] type = Diffusion variable = pressure_F [../]
[]

# [Materials]
# [flowAndTransport]
# type = FlowAndTransport 
# #block = 0
# k = 1e-6
# phi = 0.2 
# kFrac = '1e-1 1e-5'
# phiFrac = 0.4
# pressure = pressure
# inclusions_list = inclusionsList # inclusionsListRefinement
# []
# []

[BCs]
[dirinflow] type = DirichletBC variable = pressure_F value = 4 boundary = 1 []
[diroutflow] type = DirichletBC variable = pressure_F value = 1 boundary = 2 []
[]


[Preconditioning]
[./SMP]
type = SMP
full = true
[../]
[]

[Executioner]

type = Transient
solve_type= LINEAR
line_search = none
petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
petsc_options_value='  preonly   lu       NONZERO               mumps '
petsc_options = '-ksp_monitor_singular_value '
dt = 1.0
start_time = 0.0
end_time = 1.0

 #[Quadrature] type = GRID   order = ELEVENTH []

[]

[AuxKernels]
  [./initial_cond_aux]
    type = SolutionAux
    solution = soln
    execute_on = initial
    variable = aux_sol_F
  [../]
[]


[Outputs]
file_base = AdvectionOut_flow_frac
exodus = true
[]



[UserObjects]
# [./my]
# type = SolveDiffusion
# execute_on = initial
# material_name = flowAndTransport
# dirichlet_nodeset_names  = ' 11 22 '
# dirichlet_function_names = ' gamma1      gamma2 '
# neumann_sideset_names  = ' '
# neumann_function_names = ' '
# solver_type = 3
# # nemesis_file = outputUO.e
# exodus_file  = outputUO2_400_7_n.e
# variable_name = 'pressure2'
# [../]

# [./storeOperatorsUO]
# type = StoreOperators
# [../]

# [./assFID]
# type = ComputeFluxID
# execute_on = 'timestep_begin'
# block_0_id='1'
# block_1_id='2'
# material_name = flowAndTransport
# dirichlet_nodeset_names  = '2 3'
# dirichlet_function_names = 'gamma1 gamma1'
# neumann_sideset_names  = '1'
# neumann_function_names = 'gamma2'
# variable_name_1='flux_1'
# variable_name_2='flux_2'
# material_value_p = '1 1' 
# material_block_id = '1 2'
# exodus_file=flux.e
# #mesh=porous_matrix_out.e
# [../]

# [./assF]
# type = ComputeFlux
# execute_on = 'timestep_end'
# block_0_id='1'
# block_1_id='2'
# material_name = flowAndTransport
# dirichlet_nodeset_names  = ' 11 22 '
# dirichlet_function_names = ' gamma1      gamma2 '
# neumann_sideset_names  = ' '
# neumann_function_names = ' '
# variable_name_1='flux_1'
# variable_name_2='flux_2'
# #exodus_file=flux.e
# mesh=porous_matrix_out.e
# [../]

[./soln]
    type = SolutionUserObject
    mesh = fracture_network_out.e
    timestep = LATEST
    system_variables = pressure
    execute_on = 'initial'
[../]
[]



[Functions]
[gamma1]  type = ParsedFunction value = 1.0 [../]
[gamma2]  type = ParsedFunction value = 1.0 [../]
[]