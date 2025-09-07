function build(instance::String; dir=joinpath(dirname(@__DIR__), "instances"))
    
    df = CSV.read("$dir/$instance.vrp", DataFrame, silencewarnings=true)
    
    k₁ = findfirst(contains("DIMENSION"), df[:,1])
    k₂ = findfirst(contains("CAPACITY"), df[:,1])
    k₃ = findfirst(contains("NODE_COORD_SECTION"), df[:,1])
    k₄ = findfirst(contains("DEMAND_SECTION"), df[:,1])

    # nodes
    N = Vector{Node}()
    for i in 1:(parse(Int, df[k₃, 2]) - 1)
        x = parse(Int, df[k₃ + i, 2])
        y = parse(Int, df[k₃ + i, 3])
        q = parse(Float64, df[k₄ + i, 2])
        push!(N, Node(i, x, y, q, 0, 0, 0))
    end

    # arcs
    A = Matrix{Arc}(undef, length(N), length(N))
    for i in 1:length(N)
        for j in 1:length(N)
            if i != j
                l = sqrt((N[i].x - N[j].x)^2 + (N[i].y - N[j].y)^2)
                A[i, j] = Arc(i, j, l)
            else
                A[i, j] = Arc(i, j, 0.0)
            end
        end
    end
    # vehicles
    V = Vector{Vehicle}()
    for i in 1:parse(Int, df[k₂, 2])
        q = parse(Float64, df[k₂, 2])
        push!(V, Vehicle(i, 1, 1, q, 0, 0.0))
    end
    # Graph
    G = (N, A, V)
    return G
end