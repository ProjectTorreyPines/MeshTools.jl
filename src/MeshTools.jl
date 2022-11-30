module MeshTools

using LinearAlgebra

import Gmsh
using Meshes
using FileIO, MeshIO

using Memoize, LRUCache

include("mesh.jl")
export create_mesh

include("interpolation.jl")
export interpolate, barycentric_coordinates

include("integration.jl")
export integrate, quad_points


end # module MeshTools
