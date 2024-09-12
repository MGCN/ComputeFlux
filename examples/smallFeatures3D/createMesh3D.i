[Problem]
solve = false
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
 outputFileName = mesh.e
[]


[]

[Executioner]
 type = Steady
 solve_type = LINEAR
 petsc_options_iname = '-pc_type -pc_hypre_type '
 petsc_options_value = ' hypre    boomeramg     '
[]
