[Problem]
solve = false
[]

[Mesh]
[gmg]
 type = GeneratedMeshGenerator
 dim = 2
 nx = ${nxStringMoose}
 ny = ${nxStringMoose}
[]

[inclusionsList]
input = gmg 
type = InclusionsMeshModifier
fn = 2
fx_string = '0.5 0.5'
fy_string = '0.5 0.5'
fd1_string = '0.8 0.8'
fd2_string = '0.05 0.05'
fa1_string = '0.0 30.0'
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
