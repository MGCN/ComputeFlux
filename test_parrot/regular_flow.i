[Mesh]
[gmg]
type = FileMeshGenerator
file = porous_matrix_0_out.e
[]


# [./subdomains_0]
# input = gmg
# type = SubdomainBoundingBoxGenerator
# bottom_left = '0.0   0.0  0.0'
# top_right =   '0.49999   1.0001   0.0'
# block_id = 1
# [../]


# [./subdomains_1]
# input = subdomains_0
# type = SubdomainBoundingBoxGenerator
# bottom_left = '0.49999  0.000 0.0'
# top_right =   '1.0001 1.00001 0.0'
# block_id = 2
# [../]

[./subdomains_0]
input = gmg
type = SubdomainBoundingBoxGenerator
bottom_left = '0.0   0.0  0.0'
top_right =   '1.0   0.7  0.0'
block_id = 1
[../]


[./subdomains_1]
input = subdomains_0
type = SubdomainBoundingBoxGenerator
bottom_left = '0.0   0.7 0.0'
top_right =   '1.0001 1.0001 0.0'
block_id = 2
[../]

[]

[Problem]
type = FEProblem
solve=false
[]


[Variables]
[pressure_M] 
 [../]
[]

[AuxVariables]
[flux_1]
[../]
[flux_2]
[../]
[aux_sol_M]
[../]
[lambda]
[../]
[]
 

[Kernels]
[permeabilityDiffusion] type = Diffusion variable = pressure_M [../]
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
[dirinflow] type = DirichletBC variable = pressure_M value = 4 boundary = 1 []
[diroutflow] type = DirichletBC variable = pressure_M value = 1 boundary = 2 []
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

[Quadrature] type = GRID   order = ELEVENTH []

[]

[AuxKernels]
[./initial_cond_aux]
type = SolutionAux
solution = soln
execute_on = initial
variable = aux_sol_M
[../]
[]


[Outputs]
file_base = AdvectionOut_flow
exodus = true
[]


[MultiApps]
 active='transfer_1'
[./transfer_1]
 type =  TransientMultiApp
 app_type=parrot2App
 execute_on=timestep_end
 input_files=regular_flow_fracture.i
 positions='0.0 0.0 0.0'
[../]
[]


[UserObjects]
[./assFID]
type = ComputeFluxUtopiaHD
execute_on = 'timestep_end'
block_0_id='1'
block_1_id='2'
material_name = flowAndTransport
dirichlet_nodeset_names_M  = '2'
dirichlet_function_names_M = 'gamma1'
dirichlet_nodeset_names_F  = '2 3'
dirichlet_function_names_F = 'gamma1 gamma1'
neumann_sideset_names_M  = '1'
neumann_function_names_M = 'gamma2'
variable_name_1='flux_1'
variable_name_2='flux_2'
material_value_m = '1 1' 
material_value_f = '1.0'
material_block_id = '1 2'
#exodus_file=flux.e
multi_app = transfer_1
material=false
[../]


# [UserObjects]
# [./assFID]
# type = ComputeFluxUtopiaHD
# execute_on = 'timestep_end'
# block_0_id='1'
# block_1_id='2'
# material_name = flowAndTransport
# dirichlet_nodeset_names_M  = '2'
# dirichlet_function_names_M = 'gamma1'
# dirichlet_nodeset_names_F  = '2 3'
# dirichlet_function_names_F = 'gamma2 gamma2'
# neumann_sideset_names_M  = '1'
# neumann_function_names_M = 'gamma1'
# neumann_sideset_names_F  = '1'
# neumann_function_names_F = 'gamma2'
# variable_name_1='flux_1'
# variable_name_2='flux_2'
# material_value_m = '1 1' 
# material_value_f = '1e4'
# material_block_id = '1 2'
# multi_app = transfer_1
# material=false
# [../]


[./soln]
    type = SolutionUserObject
    mesh = porous_matrix_0_out.e
    timestep = LATEST
    system_variables = pressure
    execute_on = 'initial'
[../]
[]



[Functions]
[gamma1]  type = ParsedFunction value = 1.0 [../]
[gamma2]  type = ParsedFunction value = 1.0 [../]
[gamma3]  type = ParsedFunction value = 200.0 [../]
[]