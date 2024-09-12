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
fn = 10
fx_string = ' 0.135 0.15 0.3 0.275 0.74986 0.77486 0.725 0.575 0.85 0.275'
fy_string = ' 0.2392 0.205 0.36 0.70835 0.50046 0.20131 0.32375 0.84285 0.88645 0.9045'
fd1_string = ' 0.39234 0.24413 0.61774 0.48594 0.69499 0.16418 0.27415 0.51827 0.24523 0.28479'
fd2_string = ' 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001'
fa1_string = ' 115.6768 145.008 119.0546 120.9617 106.7009 155.7724 155.7723 150.2592 144.6443 208.6169'
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
