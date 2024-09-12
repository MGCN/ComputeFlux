[Mesh]
[gmg]
 type = DistributedRectilinearMeshGenerator # GeneratedMeshGenerator #
 dim = 2
 nx = ${nxMoose}
 ny = ${nxMoose}
# elem_type = QUAD9
[]

[inclusionsList]
input = gmg 
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
 doBoundaryRefinement = true
 refinements = '${refMoose} 0'
 # outputFileName = test.e
[]

[]

[Variables]
[pressure] []
# [pressure] order = THIRD family = HERMITE []
# [pressure] order = FIFTH family = HIERARCHIC []
[]

[Kernels]
[./diffusion] type = PermeabilityDiffusion variable = pressure [../]
[]

[Materials]
[./flowAndTransport] type = FlowAndTransport k = 1
kFrac = ${permMoose}
inclusions_list = inclusionsList
[../]
[]

[BCs]
[./left] type = NeumannBC   variable = pressure value = 1.0 boundary = left  [../]
[./rite] type = DirichletBC variable = pressure value = 1.0 boundary = right [../]
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
 file_base = DiffusionOut_${caseMoose}_${nxMoose}_${refMoose} 
 nemesis = true
[]
