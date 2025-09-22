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
    for i ∈ 1:n
        x = parse(Int, split(df[k₃+i,1])[2])
        y = parse(Int, split(df[k₃+i,1])[3])
        q = parse(Int, split(df[k₄+i,1])[2])
        N[i] = Node(i, x, y, q)
    end
    # arcs
    A = Matrix{Arc}(undef, n, n)
    for i ∈ 1:n
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
    q = parse(Int, df[k₂,2])
    m = d ÷ q + 1
    V = Vector{Vehicle}(undef, m)
    for i ∈ 1:m V[i] = Vehicle(i, q) end
    # create graph
    G = Graph(N, A, V)
    return G
end

"""

"""
# TODO
function clarke_wright_savings(G)
    s = Solution(G)
    G = s.G
    N = G.N
    A = G.A
    V = G.V
    
    K = eachindex(N)
    q = V[1].q

    # Initialise each customer in its own route
    # The depot is node 1

    # TODO: If you have insertnode function, we can directly use that here instead of doing all these updates
    for i ∈ K
        if isone(i) continue end 
        n = s.N[i]
        v = s.V[i-1]
        n.v = i - 1                # Each customer assigned to a unique vehicle
        n.t = 1                    # depot is tail node
        n.h = 1                    # depot is head node
        v.l = s.N[i].q             # vehicle load is customer demand
        v.n = 1                    # one customer per vehicle
        # TODO
        #=
        v.x = ?
        v.y = ?
        s.c = ?
        =#
    end

    # Calculate savings for each pair of customers
    Δ = [] # TODO: This is type unstable vector, defined as, Δ = Any{}[]. Instead use Δ as a matrix, defined as Δ = zeros(Float64, (K,K))
    for i ∈ K
        if isone(i) continue end
        for j ∈ K
            if isone(j) continue end
            if isequal(i,j) continue end
            δ = A[i,1].c + A[1,j].c - A[i,j].c
            push!(Δ, (δ, i, j)) # TODO: Consequently, this changes to Δ[i,j] = δ
        end
    end

    # TODO: Update this for matrix sort
    # Sort savings in descending order
    Δ = sort!(Δ, by = x -> -x[1], rev = true) # descending order preseves forward stability

    # TODO: Use removennode and insertnode functions here
    # Merge routes based on savings
    for (s, i, j) in Δ
        vi = s.N[i].v
        vj = s.N[j].v

        if vi != vj && s.V[vi].l + s.V[vj].l <= q  # check if different routes and merge isf feasible
            if s.N[i].t == 1 && s.N[j].h == 1      # i is at the end of its route and j is at the start of its route

                # Merge route of j into route of i
                s.v[vi].l += s.V[vj].l
                s.V[vi].n += s.V[vj].n

                # Update 
                s.N[i].h = j
                s.N[j].t = i

                # TODO: This reassigns only one node, when in fact we need to merge the two routes
                # Reassign nodes in route j to route i
                for n in s.N
                    if n.v == vj
                        n.v = vi
                    end
                end

                # Reset vehicle vj
                s.V[vj].l = 0
                s.V[vj].n = 0
            end
        end
    end

    # TODO: If you replace the above code with removenode and insertnode, these calculalations will not be needed.
    # calculalate the total cost
    total_cost = 0.0
    for v in solution.V
        if v.l > 0        # only consider active vehicles
           path = [1]    # start from depot
           current = 1
           first_node = -1
            for k in 2:n 
                if solution.N[k].v == v.i && solution.N[k].t == current
                    first_node = k
                    break
                end
            end
            if first_node == -1
            push!(path, first_node)
            current = first_node
                while true 
                next_node = solution.N[current].h
                    if next_node == 1
                        break
                    end
                end
                push!(path, next_node)
                current = next_node
            end
        end
        push!(path, 1)  # return to depot

        # Calculate cost of arcs
        for k in 1:length(path)-1
            i = path[k]
            j = path[k+1]
            total_cost += A[i,j].c
        end
    end
    solution.c = total_cost
    return solution
end