function build(instance::String; dir=joinpath(dirname(@__DIR__), "instances"))
    # read instance file
    df = CSV.read("$dir/$instance.vrp", DataFrame, silencewarnings=true)
    # fetch key indices
    k₁ = findfirst(contains("DIMENSION"), df[:,1])
    k₂ = findfirst(contains("CAPACITY"), df[:,1])
    k₃ = findfirst(contains("NODE_COORD_SECTION"), df[:,1])
    k₄ = findfirst(contains("DEMAND_SECTION"), df[:,1])
    # nodes
    n = parse(Int, df[k₁,2])
    N = Vector{Node}(undef, n)
    # TODO: check outcomes for x, y, and q
    for i in 1:n
        x = parse(Int, df[k₃+i, 2])
        y = parse(Int, df[k₃+i, 3])
        q = parse(Float64, df[k₄+i, 2])
        N[i] = Node(i, x, y, q, 0, 0, 0)
    end
    # arcs
    A = Matrix{Arc}(undef, n, n)
    for i in 1:n
        xᵢ = N[i].x
        yᵢ = N[i].y 
        for j in 1:n
            xⱼ = N[j].x
            yⱼ = N[j].y 
            c = ((xⱼ-xᵢ)^2 + (yⱼ-yᵢ)^2)^0.5
            A[i,j] = Arc(i, j, c)
        end
    end
    # vehicles
    d = 0
    for i ∈ 1:n d += N[i].q end
    m = d ÷ q + 1
    V = Vector{Vehicle}(undef, m)
    for i in 1:m V[i] = Vehicle(i, q) end
    # create graph
    G = (N, A, V)
    return G
end