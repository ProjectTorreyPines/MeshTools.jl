module MeshTools

using LinearAlgebra

import Gmsh
using Meshes
using FileIO, MeshIO

include("mesh.jl")
export create_mesh

include("interpolation.jl")
export interpolate, barycentric_coordinates
#export integrate


end # module MeshTools
