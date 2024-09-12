[Problem]
# solve = false
# skip_nl_system_check = true
[]

[Problem]
type = ParrotProblem3
use_AFC = true
operator_userobject = storeOperatorsUO
solver_type = 1
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
fn = 10
fx_string = '0.5 0.5 0.5 0.5 0.2 0.2 0.77 0.83                0.5            0.5'
fy_string = '1.125 0.175 1.6 1.6 2.05 2.05 2.05 2.05          2.25           2.25'
fz_string = '0.5 0.5 0.675 0.31 0.5 0.5 0.5 0.5               0.333333333333 0.66666666666'
fd1_string = '0.9 0.9 0.9 0.9 0.4 0.4 0.4 0.4                 1.0            1.0'
fd2_string = '1.75 0.25 1.25 1.2472 0.30594 0.30594 0.3 0.3   0.0001         0.0001'
fd3_string = '0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01         0.0001         0.0001'
fa1_string = '0 0 0 0 78.6901 -78.6901 0 0                    0.0            0.0'
fa2_string = '0 90 0 0 -90 -90 -90 -90                        0.0            0.0'
fa3_string = '0 0 16.2602 -15.8192 90 -90 0 0                 0.0            0.0'
[]

[inclusionRefinement]
 input = inclusionsList
 type  = InclusionRefinement
 inclusions_list = inclusionsList
 doBoundaryRefinement = true
 refinements = '2 0'
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

[nodeset3]
 input = nodeset2
 type = BoundingBoxNodeSetGenerator
 new_boundary =  11
  bottom_left = '-0.00001 0.0  0.333333'
  top_right  = '1.00001  0.0  0.666667'
[]
[]



[Variables]
[CM] []
[]

[AuxVariables]
[pressure2] []
[]

[Kernels]
[upwind] type = Advection variable = CM [../]
[./time] type = PorosityTimeDerivative variable = CM lumping = true [../]
[]

[Materials]
[./flowAndTransport] type = FlowAndTransport k = 1.0 kFrac = 1.0e4 phi = 0.2 phiFrac = 0.2 
pressure=pressure2 inclusions_list = inclusionsList [../]
[]

[BCs]
[./u_injection_left] type = DirichletBC boundary = 11  variable = CM value=1 [../]
[]



[Executioner]

 type = Transient
 solve_type= LINEAR
 line_search = none
 


 dt = 0.01
 num_steps=2

 [Quadrature] type = GRID   order = SEVENTH []



[]

[UserObjects]
[./my]
 type = SolveDiffusion
 execute_on = initial
 material_name = flowAndTransport
 dirichlet_nodeset_names  = ' gamma_out_1 gamma_out_2 '
 dirichlet_function_names = ' gamma1      gamma2      '
 neumann_sideset_names  = 'bottom '
 neumann_function_names = 'bottom3 '
 solver_type = 3
 # nemesis_file = outputUO.e
 exodus_file  = outputUO.e
 variable_name = 'pressure2'
[../]

[./storeOperatorsUO]
type = StoreOperators
[../]

[./MassAssembly]
type = AssembleMassMatrix
material_name = flowAndTransport
operator_userobject = storeOperatorsUO 
execute_on = 'initial'
constrain_matrix = true
dc_boundaries ='11'
dc_variables='CM'
value_D_bc='1.0'
[../]

[./assembleVolumeVectors]
type=AssembleVolumeVectors
inclusions_list = inclusionsList
execute_on = 'initial'
[../]

[]

[Postprocessors]


[./volume0]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 0
inclusions_list = inclusionsList
#  execute_on = 'initial'
[../]


[./volume1]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 1
  inclusions_list = inclusionsList
#  execute_on = 'initial'
[../]

[./volume2]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 2
  inclusions_list = inclusionsList
#  execute_on = 'initial'
[../]


[./volume3]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 3
  inclusions_list = inclusionsList
#  execute_on = 'initial'
[../]

[./volume4]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 4
  inclusions_list = inclusionsList
#  execute_on = 'initial'
[../]


[./volume5]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 5
inclusions_list = inclusionsList
#execute_on = 'initial'
[../]


[./volume6]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 6
  inclusions_list = inclusionsList
#  execute_on = 'initial'
[../]




[./volume7]
type = ElementIntegralVolumePostprocessor
fractureRegionId = 7
inclusions_list = inclusionsList
#  execute_on = 'initial'
[../]

[./Concentration0]
type = ElementIntegralConcentrationPostprocessor
variable = CM
fractureRegionId = 0
inclusions_list = inclusionsList
[../]


[./v0fast] type = IntegralSolutionOverRegionFast region = 0 doDomainSize = 1 VolumeUserObject = assembleVolumeVectors [../]
[]


[Functions]
 [gamma1]  type = ParsedFunction value = 0.0 []
 [gamma2]  type = ParsedFunction value = 0.0 []
 [bottom3] type = ParsedFunction value = (1.0/3.0<=z)*(z<=2.0/3.0) []
 #[bottom3] type = ParsedFunction value = 1.0 []
[]

[Outputs]
 file_base = output
 exodus = true
[]

