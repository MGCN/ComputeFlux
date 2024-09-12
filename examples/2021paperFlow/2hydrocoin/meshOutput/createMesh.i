[Problem]
solve = false
[]

[Mesh]
[gmg]
 type = FileMeshGenerator
 file = hydrocoin.e
[]

[inclusionsList]
 input = gmg
 type = InclusionsMeshModifier
 fn = 2
 fx_string = '0.0 0.0'
 fy_string = '-6500.0 500.0'
 fa1_string = '79.695153531233970 -45.0'
 fd1_string = '20000000.0 200000000.0'
 fd2_string = '14.758048690536125 7.071067811865476'
[]

[inclusionRefinement]
 input = inclusionsList
 type  = InclusionRefinement
 inclusions_list = inclusionsList
 doBoundaryRefinement = true
 refinements = '${refStringMoose} 0'
 outputFileName = mesh_${refStringMoose}.e
[]


[]


[Executioner]
type = Steady
[]
