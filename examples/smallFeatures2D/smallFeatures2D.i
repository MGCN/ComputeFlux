[Problem]
# solve = false
[]

[Mesh]
 [gmg]
 type = DistributedRectilinearMeshGenerator # GeneratedMeshGenerator
 dim = 2
 nx = ${nxMoose}
 ny = ${nyMoose}
 ymax = 2.25
[]

[nodeset1]
 input = gmg 
 type = BoundingBoxNodeSetGenerator
 new_boundary = gamma_out_1
 top_right =   '0.66666 2.24999 0.0'
 bottom_left = '1.00001 2.25001 0.0'
[]

[nodeset2]
 input = nodeset1 
 type = BoundingBoxNodeSetGenerator
 new_boundary = gamma_out_2
 top_right =   ' 0.333334 2.24999 0.0'
 bottom_left = '-0.000001 2.25001 0.0'
[]

[inclusionsList]
 input = nodeset2
 type = InclusionsMeshModifier
 fn = 3
fx_string = '0.5  0.675 0.31 '
fy_string = '1.125  1.6 1.6 '
fd1_string = '0.01  0.01 0.01 '
fd2_string = '1.75  1.25 1.2472 '
fa1_string = '0  -16.2602 15.8192'
[../]

#[inclusionRefinement]
# input = inclusionsList
# type  = InclusionRefinement
# inclusions_list = inclusionsList
# doBoundaryRefinement = true
# refinements = '2 0'
# outputFileName = test.e
#[]

[]

[Variables]
[pressure] []
[]

[Kernels]
[./diffusion] type = PermeabilityDiffusion variable = pressure [../]
[]

[Materials]
[./flowAndTransport] type = FlowAndTransport k = 1.0 kFrac = 1.0e4 inclusions_list = inclusionsList [../]
[]

[BCs]
[./inflow]   type = FunctionNeumannBC variable = pressure boundary = bottom      function = (x>0.33333333)*(x<0.66666667) [../]
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
 file_base = output_${nxStringMoose} 
 nemesis = true
[]
