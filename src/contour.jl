function _get_neighboring_point(pts, face, edge, d21, c12)
    edges = filter(x -> x != edge, d21(face))

    next_edge = 0
    for ed in edges
        if ed in keys(pts)
            next_edge = ed
            break
        end
    end
    if next_edge == 0
        return NTuple{2}[], edge
    end
    x, y = coordinates(pts[next_edge])
    next_pt = [(x, y)]
    delete!(pts, next_edge)

    next_face = filter(x -> x != face, c12(next_edge))
    if length(next_face) == 0
        return next_pt, next_edge
    else
        c, last_edge = _get_neighboring_point(pts, next_face[1], next_edge, d21, c12)
        return append!(next_pt, c), last_edge
    end
end

"""
    contour(f::F, m::Mesh, l::T) where {F<:Function, T<:Real) -> Vector{Vector{NTuple{2}}}

Find contour lines of function `f` at level `l`. The function `f` can be of the form f(x,y) or f(p::Point2).
"""
function contour(f::F, m::Mesh, l::T) where {F<:Function,T<:Real}

    verts = collect(vertices(m))
    elems = collect(elements(m))

    !(eltype(elems) <: Triangle) && error("Contour only supports pure triangular meshes")

    ff(x) = applicable(f, x) ? f(x) : f(coordinates(x)...)
    fvals = [ff(p) for p in verts]
    fio = fvals .> l

    topo = convert(HalfEdgeTopology, topology(m))
    d21 = Boundary{2,1}(topo) # edges in element i
    d10 = Boundary{1,0}(topo) # points in edge i
    c12 = Coboundary{1,2}(topo) # elements that have edge i


    # calculate all intersection points on the edges
    N_edges = nfacets(topo)
    edge_pts = Dict{Int,Point2}()
    for i in 1:N_edges
        i1, i2 = d10(i) # points on edge
        if xor(fio[i1], fio[i2])
            t = (l - fvals[i1]) / (fvals[i2] - fvals[i1])
            ipt = verts[i1] + t * (verts[i2] - verts[i1])
            edge_pts[i] = ipt
        end
    end

    cs = Vector{NTuple{2}}[]
    while !isempty(edge_pts)
        first_edge, first_pt = first(edge_pts)
        first_faces = c12(first_edge)
        delete!(edge_pts, first_edge)
        fp = [(coordinates(first_pt)...,)]

        c1, last_edge = _get_neighboring_point(edge_pts, first_faces[1], first_edge, d21, c12)
        c = append!(fp, c1)

        if length(first_faces) > 1
            c2, last_edge = _get_neighboring_point(edge_pts, first_faces[2], first_edge, d21, c12)
            c = append!(reverse!(c2), c)
        end
        if first_edge == last_edge
            push!(c, first(c))
        end
        push!(cs, c)
    end

    return cs
end
