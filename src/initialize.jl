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
        x = parse(Int, split(df[k₃+i,1])[2])
        y = parse(Int, split(df[k₃+i,1])[3])
        q = parse(Int, split(df[k₄+i,1])[2])
        N[i] = Node(i, x, y, q)
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
    q = parse(Int, df[k₂,2])
    m = d ÷ q + 1
    V = Vector{Vehicle}(undef, m)
    for i in 1:m V[i] = Vehicle(i, q) end
    # create graph
    G = (N, A, V)
    return G
end

"""

"""
function clarke_wright_savings(G)
    N, A, V = G
    solution = Solution(G...)
    n = length(N)
    q = V[1].q

    # Initialise each customer in its own route
    # The depot is node 1

    for i in 2:n
        solution.N[i].v = i - 1                # Each customer assigned to a unique vehicle
        solution.N[i].t = 1                    # depot is tail node
        solution.N[i].h = 1                    # depot is head node
        solution.V[i-1].l = solution.N[i].q    # vehicle load is customer demand
        solution.V[i-1].n = 1                  # one customer per vehicle
    end

    # Calculate savings for each pair of customers
    savings = []
    for i in 2:n
        for j in i+1:n
            s = A[i,1].c + A[1,j].c - A[i,j].c
            push!(savings, (s, i, j))
        end
    end

    # Sort savings in descending order
    savings = sort!(savings, by = x -> -x[1], rev = true) # descending order preseves forward stability

    # Merge routes based on savings
    for (s, i, j) in savings
        vi = solution.N[i].v
        vj = solution.N[j].v

        if vi != vj && solution.V[vi].l + solution.V[vj].l <= q  # check if different routes and merge isf feasible
            if solution.N[i].t == 1 && solution.N[j].h == 1      # i is at the end of its route and j is at the start of its route

                # Merge route of j into route of i
                solution.v[vi].l += solution.V[vj].l
                solution.V[vi].n += solution.V[vj].n

                # Update 
                solution.N[i].h = j
                solution.N[j].t = i

                # Reassign nodes in route j to route i
                for node in solution.N
                    if node.v == vj
                        node.v = vi
                    end
                end

                # Reset vehicle vj
                solution.V[vj].l = 0
                solution.V[vj].n = 0
            end
        end
    end

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