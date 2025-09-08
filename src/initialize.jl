function build(instance::String; dir=joinpath(dirname(@__DIR__), "instances"))
    # read instance file
    df = CSV.read("$dir/$instance.vrp", DataFrame, silencewarnings=true)
    # fetch key indices
    k₁ = findfirst(contains("DIMENSION"), df[:,1])
    k₂ = findfirst(contains("CAPACITY"), df[:,1])
    k₃ = findfirst(contains("NODE_COORD_SECTION"), df[:,1])
    k₄ = findfirst(contains("DEMAND_SECTION"), df[:,1])
    # nodes
    # NOTE: made some improvements (remove this comment once you are satisfied with the changes)
    n = parse(Int, df[k₁,2])
    N = Vector{Node}(undef, n)
    for i in 1:n
        x = parse(Int, df[k₃+i, 2])
        y = parse(Int, df[k₃+i, 3])
        q = parse(Float64, df[k₄+i, 2])
        N[i] = Node(i, x, y, q, 0, 0, 0)
    end
    # arcs
    # NOTE: made some improvements (remove this comment once you are satisfied with the changes)
    A = Matrix{Arc}(undef, n, n)
    for i in 1:n
        xᵢ = N[i].x
        yᵢ = N[i].y 
        for j in 1:n
            xⱼ = N[j].x
            yⱼ = N[j].y 
            l  = ((xⱼ-xᵢ)^2 + (yⱼ-yᵢ)^2)^0.5
            A[i,j] = Arc(i, j, l)
        end
    end
    # vehicles
    # TODO: This part is incorrect and inefficient
    V = Vector{Vehicle}()
    for i in 1:parse(Int, df[k₂, 2])
        q = parse(Float64, df[k₂, 2])
        push!(V, Vehicle(i, 1, 1, q, 0, 0.0))
    end
    # create graph
    G = (N, A, V)
    return G
end