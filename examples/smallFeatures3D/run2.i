[Problem]
# solve = false
# skip_nl_system_check = true
[]

[Mesh]
[gmg]
 type = GeneratedMeshGenerator  # DistributedRectilinearMeshGenerator #
 dim = 3
 xmin = 0.0
 xmax = 1.0
 ymin = 0.0
 ymax = 2.25
 zmin = 0.0
 zmax = 1.0
 nx = 12
 ny = 27
 nz = 12
[]

[inclusionsList]
input = gmg
type = InclusionsMeshModifier
fn = 8
fx_string = '0.5 0.5 0.5 0.5 0.2 0.2 0.77 0.83'
fy_string = '1.125 0.175 1.6 1.6 2.05 2.05 2.05 2.05'
fz_string = '0.5 0.5 0.675 0.31 0.5 0.5 0.5 0.5'
fd1_string = '0.9 0.9 0.9 0.9 0.4 0.4 0.4 0.4'
fd2_string = '1.75 0.25 1.25 1.2472 0.30594 0.30594 0.3 0.3'
fd3_string = '0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01'
fa1_string = '0 0 0 0 78.6901 -78.6901 0 0'
fa2_string = '0 90 0 0 -90 -90 -90 -90'
fa3_string = '0 0 16.2602 -15.8192 90 -90 0 0'
[]

[inclusionRefinement]
 input = inclusionsList
 type  = InclusionRefinement
 inclusions_list = inclusionsList
 doBoundaryRefinement = true
 refinements = '5 0'
# outputFileName = mesh.e
[]

[nodeset1]
 input = inclusionRefinement 
 type = BoundingBoxNodeSetGenerator
 new_boundary = gamma_out_1
 top_right =   '-0.00001 2.24999 -0.00001'
 bottom_left = ' 1.00001 2.25001  0.33334'
[]

[nodeset2]
 input = nodeset1 
 type = BoundingBoxNodeSetGenerator
 new_boundary = gamma_out_2
 top_right =   '-0.00001 2.24999  0.66666'
 bottom_left = ' 1.00001 2.25001  1.00001'
[]

[]

[Variables]
[pressure] []
[]

[Kernels]
[./diffusion] type = PermeabilityDiffusion variable = pressure [../]
[]

[Materials]
[./flowAndTransport] type = FlowAndTransport k = 1.0 kFrac = 1.0e2 inclusions_list = inclusionsList [../]
[]

[BCs]
[./inflow]   type = FunctionNeumannBC variable = pressure boundary = bottom      function = (z>0.33333333)*(z<0.66666667) [../]
[./outflow1] type = DirichletBC       variable = pressure boundary = gamma_out_1 value =  0.0 [../]
[./outflow2] type = DirichletBC       variable = pressure boundary = gamma_out_2 value =  0.0 [../]
[]

[Executioner]
 type = Steady
 solve_type = LINEAR
 petsc_options_iname = '-pc_type -pc_hypre_type '
 petsc_options_value = ' hypre    boomeramg     '
[]

[Outputs]
 file_base = output
 # nemesis = true
 exodus = true
[]
