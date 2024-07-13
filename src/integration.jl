"""
    quad_points(tri::Triangle, ::Val{N}; w=1.0) where N

Given a triangle
returns a list of evaluation points [(r,z,weight),...]
n   number of quadrature points.
currently: 1, 3 or 6
w   weight scale factor default = 1
Coefficients taken from http://www.cs.rpi.edu/~flaherje/pdf/fea6.pdf
Joseph E. Flaherty course notes, Rensselaer Polytechnic Institute
"""
function quad_points(tri::Triangle, ::Val{N}; w=1.0) where {N}

    w *= measure(tri)
    verts = vertices(tri)
    r1, z1 = coordinates(verts[1])
    r2, z2 = coordinates(verts[2])
    r3, z3 = coordinates(verts[3])

    if N == 1
        # One point in the middle of the triangle
        return [((r1 + r2 + r3) / 3, (z1 + z2 + z3) / 3, 1.0 * w)]
    elseif N == 3
        return [
            ((4 * r1 + r2 + r3) / 6, (4 * z1 + z2 + z3) / 6, w * 1.0 / 3),
            ((r1 + 4 * r2 + r3) / 6, (z1 + 4 * z2 + z3) / 6, w * 1.0 / 3),
            ((r1 + r2 + 4 * r3) / 6, (z1 + z2 + 4 * z3) / 6, w * 1.0 / 3)
        ]
    elseif N == 6
        a = 0.816847572980459
        b = 0.5 * (1.0 - a)

        c = 0.108103018168070
        d = 0.5 * (1.0 - c)

        return [
            ((a * r1 + b * r2 + b * r3), (a * z1 + b * z2 + b * z3), w * 0.109951743655322),
            ((b * r1 + a * r2 + b * r3), (b * z1 + a * z2 + b * z3), w * 0.109951743655322),
            ((b * r1 + b * r2 + a * r3), (b * z1 + b * z2 + a * z3), w * 0.109951743655322),
            ((c * r1 + d * r2 + d * r3), (c * z1 + d * z2 + d * z3), w * 0.223381589678011),
            ((d * r1 + c * r2 + d * r3), (d * z1 + c * z2 + d * z3), w * 0.223381589678011),
            ((d * r1 + d * r2 + c * r3), (d * z1 + d * z2 + c * z3), w * 0.223381589678011)
        ]
    else
        error("Quadrature Rule not available for n=$n")
    end
end

"""
    quad_points(poly::Ngon, ::Val{N}) where N

Given a polygon
calculates a set of quadrature points and weights, by splitting
the polygon into triangles.
returns a list of evaluation points and weights [(r,z,weight),...]
These are normalised to calculate the average value of a function
over the polygon; multiply by the area to get the integral.
n   number of quadrature points in each triangle
currently: 1, 3 or 6
"""
function quad_points(poly::Ngon, ::Val{N}) where {N}

    # Split polygon into triangles
    trimesh = discretize(poly, FanTriangulation())

    quads = NTuple{3,Float64}[]  # List of all points
    for (i, t) in enumerate(elements(trimesh))
        append!(quads, quad_points(t, Val(N)))  # Quadrature points for this triangle
    end
    return quads
end

"""
    quad_points(m::Mesh, ::Val{N}) where N

Calculate quadrature points for all elements in the mesh with n order quadrature points
"""
@memoize LRU(maxsize=5) function quad_points(m::Mesh, ::Val{N}) where {N}
    quads = NTuple{3,Float64}[]
    for (i, t) in enumerate(elements(m))
        append!(quads, quad_points(t, Val(N)))
    end
    return quads
end

"""
    integrate(func, poly::T; n=6) where T<:Ngon

Integrate func(r,z) over polygon with n order quadrature points
"""
function integrate(func, poly::T; n=6) where {T<:Ngon}
    quad = quad_points(poly, Val(n))
    return sum(func(r, z) * w for (r, z, w) in quad)
end

"""
    integrate(func, m::Mesh; n=6)

Integrate func(r,z) over mesh with n order quadrature points
"""
function integrate(func, m::Mesh; n=6)
    quads = quad_points(m, Val(n))
    return sum(func(r, z) * w for (r, z, w) in quads)
end
