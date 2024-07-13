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

include("contour.jl")
export contour

const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__, all=false, imported=false) if name != Symbol(@__MODULE__)]

end # module MeshTools
