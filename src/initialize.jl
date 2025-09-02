function build(instance::String; dir=joinpath(dirname(@__DIR__), "instances"))
    
    df = CSV.read("$dir/$instance.vrp", DataFrame, silencewarnings=true)
    
    k₁ = findfirst(contains("DIMENSION"), df[:,1])
    k₂ = findfirst(contains("CAPACITY"), df[:,1])
    k₃ = findfirst(contains("NODE_COORD_SECTION"), df[:,1])
    k₄ = findfirst(contains("DEMAND_SECTION"), df[:,1])

    # Nodes

    # Arcs

    # Vehicles

    
    # Graph
    G = (N, A, V)
    return G
end