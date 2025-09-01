
function build(instance::String; dir=joinpath(@__DIR__, "..", "instances"))

    #Arc
    df = CSV.read(joinpath(dir, instance, "arcs.csv"), DataFrame, header=false)
    A = Dict{Tuple{Int,Int},Arc}()
    n = lastindex(C)
    for iᵗ in 1:n
        for iʰ in 1:n
            l = df[iᵗ, iʰ]
            a = Arc(iᵗ, iʰ, l)
            A[(iᵗ, iʰ)] = a
        end
    end

    #vehicles
    df = CSV.read(joinpath(dir, instance, "vehicles.csv"), DataFrame, header=false)
    for i in 1:size(df, 1)
        qᵛ = df[i, 1]
        π = df[i, 2]
        v = Vehicle(i, qᵛ, π, Vector{Route}())
        x = 0
        y = 0
        qᶜ= 0
        l = 0
        v.r = Route(iᵛ, x, y, qᶜ, l)
        push!(d.V, v)
    end

    #Depots
    df = CSV.read(joinpath(dir, instance, "depots.csv"), DataFrame, header=false)
    D = Vector{Depot}()
    for i in 1:size(df, 1)
        x = df[i, 1]
        y = df[i, 2]
        d = Depot(x, y, deepcopy(V))
        push!(D, d)
    end

    #customer nodes
    df = CSV.read(joinpath(dir, instance, "customers.csv"), DataFrame, header=false)
    C = Vector{Customer}()
    for i in 1:size(df, 1)
        x = df[i, 1]
        y = df[i, 2]
        qᶜ = df[i, 3]
        iᵗ = 0
        iʰ = 0
        iʳ = 0
        iᵛ = 0
        c = Customer(i, x, y, qᶜ, iᵗ, iʰ, iʳ, iᵛ)
        push!(C, c)
    end
    G = (D, C, A)
    return G
end

    
