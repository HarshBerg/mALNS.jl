"""
    visualize(instance; backend=gr)

    Plots instance using Plots.jl as frontend and GR as backend
"""
using Plots
function visualize(G; backend=gr;
                   show_routes::Union{nothing, vector{Vector{Int}}}=nothing,
                   depot_shape=:star5, depot_color=:red, depot_size=10,
                   customer_shape=:circle, customer_color=:blue, customer_size=5,
                   title::String="VRP Instance Visualization",
                   xlabel::String="X", ylabel::String="Y",
                   legend::Bool=false)

    N, A, V = G

    xs = [node.x for node in N] # or xs = [getfield(node, :x) for node in N]
    ys = [node.y for node in N] # or ys = [getfield(node, :y) for node in N]

    depot_x, depot_y = xs[1], ys[1]
    customer_xs = xs[2:end]
    customer_ys = ys[2:end]
    p = scatter(customer_x, customer_y, marker=(customer_shape, customer_size, customer_color), label="Customers")
    scatter!(p, [depot_x], [depot_y], marker=(depot_shape, depot_size, depot_color), label="Depot")

    if show_routes !== nothing
        colors = distinguishable_colors(length(show_routes))
        for (idx, route) in enumerate(show_routes)
            route_xs = [N[i].x for i in route]
            route_ys = [N[i].y for i in route]
            plot!(p, route_xs, route_ys, lw=2, line=(:solid, colors[idx]), label="Route $idx")
        end
    end
    title!(p, title)
    xlabel!(p, xlabel)
    ylabel!(p, ylabel)

    return p
end
# TODO: This is way too complex; can be significantly simplified. It also looks like a first draft of a LLM. Let's discuss it in the next meeting.

# also can add something to show lables and annotate the node indices for node in N









        