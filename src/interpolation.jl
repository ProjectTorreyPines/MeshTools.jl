function barycentric_coordinates(t::Triangle, p::Point)
    x, y = coordinates(p)
    (x1,y1),(x2,y2),(x3,y3) = coordinates.(vertices(t))

    detT = (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3)
    λ₁ =  ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3))/detT
    λ₂ =  ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3))/detT
    λ₃ = 1 - λ₁ - λ₂

    return λ₁, λ₂, λ₃
end

function interpolate(t::Triangle, F::Vector, p)
    λ₁, λ₂, λ₃ = barycentric_coordinates(t,p)

    return λ₁*F[1] + λ₂*F[2] + λ₃*F[3]
end

function interpolate(m::Mesh, vertexF::Matrix, p)
    @assert size(vals)[1] == 3
    triangles = collect(elements(m))

    ind = findfirst(t -> p ∈ t, triangles)
    if ind == nothing
        throw(DomainError(p,"Point outside of mesh"))
    end

    return interpolate(triangles[ind],vertexF[:,ind],p)
end

function interpolate(m::Mesh, vertexF::Vector, p)
    triangles = collect(elements(m))

    ind = findfirst(t -> p ∈ t, triangles)
    if ind == nothing
        throw(DomainError(p,"Point outside of mesh"))
    end

    topo = convert(HalfEdgeTopology,topology(m))
    δ₂₀ = Boundary{2,0}(topo)
    return interpolate(triangles[ind],vertexF[δ₂₀(ind)],p)
end

function interpolate(md::MeshData, key::Symbol, p)
    F = getproperty(values(md,0),key)
    m = domain(md)
    return interpolate(m, F, p)
end
