function visualize(G; backend=gr)
    backend()
    N,_ = build(instance)
    fig = plot(legend=:none)
    K = lastindex(N)
    X = zeros(Float64, K)       # abcissa
    Y = zeros(Float64, K)       # ordinate
    C = fill("color", K)        # color
    S = zeros(Int, K)           # size
    M = fill(:shape, K)         # marker
    for (k,n) ∈ enumerate(N)
        X[k] = n.x
        Y[k] = n.y
        C[k] = isdepot(n) ? "#b4464b" : "#d1e0ec"
        S[k] = isdepot(n) ? 6 : 5
        M[k] = isdepot(n) ? :rect : :circle
    end
    scatter!(X, Y, color=C, markersize=S, markershape=M, markerstrokewidth=0)
    return fig
end

function visualize(s::Solution; backend=gr)
    backend()
    # TODO
    # figure
    return fig
end







        