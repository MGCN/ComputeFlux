[Problem]
solve = false
[]

[Mesh]
[gmg]
 type = DistributedRectilinearMeshGenerator # GeneratedMeshGenerator #
 dim = 2
 nx = ${nxStringMoose}
 ny = ${nxStringMoose}
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
 refinements = '${refStringMoose} 0'
 outputFileName = mesh_${nxStringMoose}_${refStringMoose}.e
[]
[]


[Executioner]
type = Steady
[]
