"""
    create_mesh(outer::P; inner=Vector{P}[], outer_length=nothing, inner_lengths=nothing, spline=true, spline_inner=fill(spline,length(inner)) where P<:Vector{NTuple{2}}

Create a Triangular mesh from a vector of tuples using Gmsh. Interior boundaries define holes. Returns Meshes.SimpleMesh
"""
function create_mesh(outer::Vector{NTuple{2,T}}; inner=Vector{Vector{NTuple{2,T}}}[], outer_length=nothing, inner_lengths=nothing, spline=false, spline_inners=tuple(x->spline,length(inner))) where T
    if first(outer) != last(outer)
        push!(outer,first(outer))
    end

    for i=1:length(inner)
        if first(inner[i]) != last(inner[i])
            push!(inner[i],first(inner[i]))
        end
    end
    poly = PolyArea(outer, inner, fix=true)

    if outer_length == nothing
        ol = length(poly.outer)/50
    else
        ol = outer_length
    end

    if inner_lengths == nothing
        ils = ntuple(x->ol,length(inner))
    else
        ils = inner_lengths
    end

    return create_mesh(poly, ol, ils, spline, spline_inners)
end

function create_mesh(poly::PolyArea, outer_length, inner_lengths::Tuple, spline::Bool, spline_inners::Tuple)

    debug = false
    ext = poly.outer
    interior = poly.inners

    Gmsh.initialize(["-v","0"])

    name = tempname()
    Gmsh.gmsh.model.add(name)

    factory = Gmsh.gmsh.model.geo

    p_t = 0  # point tag
    l_t = 0  # line/spline/curve tag
    c_t = 0  # curve loop tag
    s_t = 0  # surface tag

    # Add interior curves
    p_t_prev = 0
    l_t_prev = 0
    for (j,h) in enumerate(interior)
        v = vertices(h)
        hN = length(v)
        for i=1:hN
            p = coordinates(v[i])
            p_t += 1
            factory.addPoint(p[1], p[2], 0.0, inner_lengths[j], p_t)
            debug && println("addPoint $(p[1]), $(p[2]), 0.0, $(inner_lengths[j]), $p_t")
        end
        if spline_inners[j]
            l_t += 1
            factory.addSpline(append!(collect((p_t_prev + 1):p_t),p_t_prev+1),l_t)
            c_t += 1
            factory.addCurveLoop([l_t],c_t)
        else
            for i=1:hN
                l_t += 1
                factory.addLine(p_t_prev + i, p_t_prev + mod(i,hN) + 1, l_t)
                debug && println("addLine $(p_t_prev+i), $(p_t_prev + mod(i,hN)+1), $l_t")
            end
            c_t += 1
            factory.addCurveLoop(collect((l_t_prev+1):l_t),c_t)
        end
        p_t_prev = p_t
        l_t_prev = l_t
    end

    # Add Exterior Surface
    v = vertices(ext)
    N = length(v)
    for i = 1:N
        p = coordinates(v[i])
        p_t += 1
        factory.addPoint(p[1], p[2], 0.0, outer_length, p_t)
        debug && println("addPoint $(p[1]), $(p[2]), 0.0, $(outer_length), $p_t")
    end
    if spline
        l_t += 1
        factory.addSpline(append!(collect((p_t_prev + 1):p_t),p_t_prev+1),l_t)
        c_t += 1
        factory.addCurveLoop([l_t],c_t)
    else
        for i=1:N
            l_t += 1
            factory.addLine(p_t_prev + i, p_t_prev + mod(i,N) + 1, l_t)
            debug && println("addLine $(p_t_prev+i), $(p_t_prev + mod(i,N)+1), $l_t")
        end
        c_t += 1
        factory.addCurveLoop(collect((l_t_prev+1):l_t),c_t)
    end
    s_t += 1 

    factory.addPlaneSurface(collect(c_t:-1:1),s_t)

    factory.synchronize()

    Gmsh.gmsh.model.mesh.generate(2)
    
    filename = name*".msh"
    Gmsh.gmsh.write(filename)
    Gmsh.finalize()
    
    m = load(filename)

    # Convert to Meshes.jl format
    points = [Tuple([p[1],p[2]]) for p in Set(m.position)]
    indices = Dict(p=>i for (i,p) in enumerate(points))
    connect = map(m) do el
        Meshes.connect(Tuple(indices[Tuple([p[1],p[2]])] for p in el))
    end
    return Meshes.SimpleMesh(points,connect)
end
