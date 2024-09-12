[Problem]
solve = false
# skip_nl_system_check = true
[]

[Mesh]
[gmg]
 type = GeneratedMeshGenerator  # DistributedRectilinearMeshGenerator #
 dim = 3
 nx = 20
 ny = 20
 nz = 20
[]

[inclusionsList]
input = gmg
type = InclusionsMeshModifier
fn = 1
fx_string = '0.5'
fy_string = '0.5'
fz_string = '0.5'
fd1_string = '0.5'
fd2_string = '0.5'
fd3_string = '0.5'
fa1_string = '0.0'
fa2_string = '0.0'
fa3_string = '0.0'
[]

[inclusionRefinement]
 input = inclusionsList
 type  = InclusionRefinement
 inclusions_list = inclusionsList
 doBoundaryRefinement = true
 refinements = '1 0'
 outputFileName = mesh.e
[]


[]

[Executioner]
 type = Steady
 solve_type = LINEAR
 petsc_options_iname = '-pc_type -pc_hypre_type '
 petsc_options_value = ' hypre    boomeramg     '
[]
