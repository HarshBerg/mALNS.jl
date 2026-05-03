"""
    build(instance::String; dir=joinpath(dirname(@__DIR__), "instances"))

Returns `Graph` for `instance` stored at `dir`.
"""
function build(instance::String; dir=joinpath(dirname(@__DIR__), "instances"))
    # read instance file
    df = CSV.read("$dir/$instance.vrp", DataFrame, silencewarnings=true)
    # fetch key indices
    k₁ = findfirst(contains("DIMENSION"), df[:,1])::Int
    k₂ = findfirst(contains("CAPACITY"), df[:,1])::Int
    k₃ = findfirst(contains("NODE_COORD_SECTION"), df[:,1])::Int
    k₄ = findfirst(contains("DEMAND_SECTION"), df[:,1])::Int
    # fetch nodes
    n = parse(Int, df[k₁,2])
    N = Vector{Node}(undef, n)
    for i ∈ 1:n
        v = ""
        try v = string(split(df[k₃+i,1])[2])
        catch e
            v = string(df[k₃+i,2])
        end
        x = parse(Int, v)
        try v = string(split(df[k₃+i,1])[3])
        catch e
            v = string(df[k₃+i,3])
        end
        y = parse(Int, v)
        try v = string(split(df[k₄+i,1])[2])
        catch e
            v = string(df[k₄+i,2])
        end
        r = isone(i) ? 0. : hypot((x - N[1].x), (y - N[1].y))
        θ = isone(i) ? 0. : atan((y - N[1].y) / (x - N[1].x))
        q = parse(Int, v)
        N[i] = Node(i, x, y, r, θ, q)
    end
    # create arcs
    A = Matrix{Arc}(undef, n, n)
    for t ∈ 1:n
        xₜ = N[t].x
        yₜ = N[t].y 
        for h ∈ 1:n
            xₕ = N[h].x
            yₕ = N[h].y
            c  = round(Int, hypot(xₕ - xₜ, yₕ - yₜ))
            a  = Arc(t, h, c)
            A[t,h] = a
        end
    end   
    # fetch vehicles
    q = parse(Int, df[k₂,2])
    V = Vector{Vehicle}(undef, n-1)
    for i ∈ 1:(n-1) V[i] = Vehicle(i, q) end
    # create graph
    G = Graph(N, A, V)
    return G
end

"""
    initialize(instance::String; dir=joinpath(dirname(@__DIR__), "instances"))

Returns initial solution for `instance` stored at `dir` using Clarke & Wright Algorithm.
"""
function initialize(instance::String; dir=joinpath(dirname(@__DIR__), "instances"))
    # pre-initialize
    G = build(instance; dir=dir)
    s = Solution(G, 0, 0)
    G = s.G
    N = G.N
    A = G.A
    V = G.V
    Ψ = (intramove!, intermove!, intraswap!, interswap!, intraopt!, interopt!)   
    L = eachindex(Ψ)
    # initialize
    K = eachindex(N)
    Δ = zeros(Int, (K,K))       # Δ[i,j]: Savings from merging node N[i] and N[j]   
    d = N[1]
    for k ∈ K
        n = N[k]
        if isdepot(n) continue end
        v = V[k-1]
        insertnode!(n, d, d, v, s)
    end
    # iterate through each node pair
    for i ∈ K
        nᵢ = N[i]
        if isdepot(nᵢ) continue end
        for j ∈ K
            nⱼ = N[j]
            if isdepot(nⱼ) continue end
            if j ≥ i continue end
            δ = (A[i,1].c + A[1,j].c) - A[i,j].c
            Δ[i,j] = δ
        end
    end
    # perform feasible greedy merge
    P = sortperm(vec(Δ), rev=true)
    T = Tuple.(CartesianIndices(Δ)[P])
    for (i,j) ∈ T
        if iszero(C[i,j]) break end
        nᵢ = N[i]
        nⱼ = N[j]
        # nodal feasibility check
        tᵢ = N[nᵢ.t]
        tⱼ = N[nⱼ.t]
        hᵢ = N[nᵢ.h]
        hⱼ = N[nⱼ.h]
        if iscustomer(tᵢ) && iscustomer(hᵢ) continue end
        if iscustomer(tⱼ) && iscustomer(hⱼ) continue end
        # vehicular feasibility check
        vᵢ = V[nᵢ.v]
        vⱼ = V[nⱼ.v]
        if isequal(vᵢ, vⱼ) continue end
        if vᵢ.l + vⱼ.l > vⱼ.q continue end
        # merge
        k = vᵢ.n
        φ = isdepot(hⱼ)
        for _ ∈ 1:k
            nᵢ = N[vᵢ.s]
            nⱼ = N[j]
            tᵢ = N[nᵢ.t]
            tⱼ = N[nⱼ.t]
            hᵢ = N[nᵢ.h]
            hⱼ = N[nⱼ.h]
            removenode!(nᵢ, tᵢ, hᵢ, vᵢ, s)
            if φ insertnode!(nᵢ, nⱼ, hⱼ, vⱼ, s)
            else insertnode!(nᵢ, tⱼ, nⱼ, vⱼ, s) end
        end
    end
    # remove non-operational vehicles
    filter!(isopt, V)
    # reset indices
    K = eachindex(V)
    for k ∈ K
        v = V[k]
        v.i = k
        n = N[v.s]
        for _ ∈ 1:v.n
            n.v = k
            n = N[n.h]
        end
    end
    # add slack vehicles
    i = lastindex(V)
    v = V[i]
    q = v.q
    for j ∈ 1:3 push!(V, Vehicle(i+j, q)) end
    # local search
    for l ∈ L localsearch!(MersenneTwister(length(N)), 200 * length(N), s, Ψ[l]) end
    # return solution
    return s
end