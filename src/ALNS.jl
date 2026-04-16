"""
    cALNS([rng::AbstractRNG], χ::ALNSparameters, sₒ::Solution; mute=false))

conventional Adaptive Large Neighborhood Search (cALNS)

Given ALNS optimization parameters `χ` and an initial solution `sₒ`, 
ALNS adaptively searches large neighborhoods in the solution domain and
returns the best found solution. Additionally, displays a convergence plot.

Takes `mute` as argument. If `true` mutes progressbar and pltcnv output.

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
@inline function removalcount(rng::AbstractRNG, s::Solution, μ̲::Float64, e̲::Int, μ̅::Float64, e̅::Int)
    customers = count(n -> iscustomer(n) && isclose(n), s.G.N)
    customers == 0 && return 0
    lo = min(customers, max(e̲, ceil(Int, μ̲ * customers)))
    hi = min(customers, max(lo, min(e̅, floor(Int, μ̅ * customers))))
    return rand(rng, lo:hi)
end

function conALNS(rng::AbstractRNG, χ::ALNSparameters, sₒ::Solution; mute=false)
    # Step 0: Pre-initialize
    j,k = χ.j, χ.k
    n,m = χ.n, χ.m
    ω̅, τ̅ = χ.ω̅, χ.τ̅
    ω̲ ,τ̲ = χ.ω̲, χ.τ̲
    μ̲ ,e̲ = χ.μ̲, χ.e̲
    μ̅ ,e̅ = χ.μ̅, χ.e̅
    θ, ρ = χ.θ, χ.ρ
    Ψᵣ, Ψᵢ, Ψₗ = χ.Ψᵣ, χ.Ψᵢ, χ.Ψₗ
    σ₁, σ₂, σ₃ = χ.σ₁, χ.σ₂, χ.σ₃
    R = eachindex(Ψᵣ)
    I = eachindex(Ψᵢ)
    L = eachindex(Ψₗ)
    X = fill(zero(UInt), 1 + j * n)
    Z = Vector{Float64}(undef, 1 + j * n)
    # Step 1: Initialize
    s = deepcopy(sₒ)
    s_best = deepcopy(sₒ)
    x = h(s)
    z = f(s)
    X[1] = x
    Z[1] = z
    z_best = z
    t = ω̅ * f(sₒ) / log(1 / τ̅)
    Cᵣ = zeros(Int, R)
    Cᵢ = zeros(Int, I)
    Pᵣ = zeros(Float64, R)
    Pᵢ = zeros(Float64, I)
    Wᵣ = ones(Float64, R)
    Wᵢ = ones(Float64, I)
    Sᵣ = zeros(Float64, R)
    Sᵢ = zeros(Float64, I)
    # Step 2: Loop over segments
    for segment ∈ 1:j
        # Step 2.1: Reset count and score for every removal and insertion operator
        for r ∈ R Cᵣ[r], Sᵣ[r] = 0, 0 end
        for i ∈ I Cᵢ[i], Sᵢ[i] = 0, 0 end
        # Step 2.2: Update selection probability for every removal and insertion operator
        for r ∈ R Pᵣ[r] = Wᵣ[r] / sum(Wᵣ) end
        for i ∈ I Pᵢ[i] = Wᵢ[i] / sum(Wᵢ) end
        # Step 2.3: Loop over iterations within the segment
        for iteration ∈ 1:n
            # Step 2.3.1: Randomly select a removal and an insertion operator based on operator selection probabilities, and consequently update count for the selected operators
            r = sample(rng, R, Weights(Pᵣ))
            i = sample(rng, I, Weights(Pᵢ))
            Cᵣ[r] += 1
            Cᵢ[i] += 1
            # Step 2.3.2: Using the selected removal and insertion operators destroy and repair the current solution to develop a new solution
            q = removalcount(rng, s, μ̲, e̲, μ̅, e̅)
            s′ = deepcopy(s)
            s′ = Ψᵣ[r](rng, q, s′)
            s′ = Ψᵢ[i](rng, s′)
            x′ = h(s′)
            z′ = f(s′)
            # Step 2.3.3: If this new solution is better than the best solution, then set the best solution and the current solution to the new solution, and accordingly update scores of the selected removal and insertion operators by σ₁
            if z′ < z_best
                s_best = s′
                s = s′
                x = x′
                z = z′
                z_best = z′
                Sᵣ[r] += σ₁
                Sᵢ[i] += σ₁
            # Step 2.3.4: Else if this new solution is only better than the current solution, then set the current solution to the new solution and accordingly update scores of the selected removal and insertion operators by σ₂
            elseif z′ < z
                s = s′
                x = x′
                z = z′
                if x′ ∉ X
                    Sᵣ[r] += σ₂
                    Sᵢ[i] += σ₂
                end
            # Step 2.3.5: Else accept the new solution with simulated annealing acceptance criterion. Further, if the new solution is also newly found then update operator scores by σ₃
            else 
                if rand(rng) < exp((z - z′) / t)
                    s = s′
                    x = x′
                    z = z′
                    if x′ ∉ X
                        Sᵣ[r] += σ₃
                        Sᵢ[i] += σ₃
                    end
                end
            end
            # Step 2.3.6: Update temperature according to cooling schedule
            t = max(ω̲  * f(s_best) / log(1 / τ̲), t * θ) 
            # Step 2.3.7: Store best solution from each iteration
            X[1 + (segment - 1) * n + iteration] = x
            Z[1 + (segment - 1) * n + iteration] = z
        end
        # Step 2.4: Update weights for every removal and insertion operator (module)
        for r ∈ R if !iszero(Cᵣ[r]) Wᵣ[r] = (1 - ρ) * Wᵣ[r] + ρ * Sᵣ[r] / Cᵣ[r] end end
        for i ∈ I if !iszero(Cᵢ[i]) Wᵢ[i] = (1 - ρ) * Wᵢ[i] + ρ * Sᵢ[i] / Cᵢ[i] end end
        # Step 2.5: Reset current solution
        s = sₒ
        # Step 2.6: Perform local search
        s′ = deepcopy(s)
        for l ∈ L
            s′ = Ψₗ[l](rng, m, s′)
        end
        x′ = h(s′)
        z′ = f(s′)
        if z′ < z_best
            s_best = s′
            s = s′
            x = x′
            z = z′
            z_best = z′
        elseif z′ < z
            s = s′
            x = x′
            z = z′
        end
    end
    # Step 3: Display the convergence plot and return the best solution
    return s_best
end



"""
    mALNS([rng::AbstractRNG], χ::ALNSparameters, sₒ::Solution; mute=false))

modular Adaptive Large Neighborhood Search (mALNS)

Given ALNS optimization parameters `χ` and an initial solution `sₒ`, 
ALNS adaptively searches large neighborhoods in the solution domain and
returns the best found solution. Additionally, displays a convergence plot.

Takes `mute` as argument. If `true` mutes progressbar and pltcnv output.

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
function modALNS(rng::AbstractRNG, χ::ALNSparameters, sₒ::Solution; mute=false)
    # Step 0: Pre-initialize
    j,k = χ.j, χ.k
    n,m = χ.n, χ.m
    ω̅, τ̅ = χ.ω̅, χ.τ̅
    ω̲ ,τ̲ = χ.ω̲, χ.τ̲
    μ̲ ,e̲ = χ.μ̲, χ.e̲
    μ̅ ,e̅ = χ.μ̅, χ.e̅
    θ, ρ = χ.θ, χ.ρ
    Ψᵣ, Ψᵢ, Ψₗ = χ.Ψᵣ, χ.Ψᵢ, χ.Ψₗ
    σ₁, σ₂, σ₃ = χ.σ₁, χ.σ₂, χ.σ₃
    Rᵣ,Rᵪ = size(Ψᵣ)
    Iᵣ,Iᵪ = size(Ψᵢ)  
    L = eachindex(Ψₗ)
    X = fill(zero(UInt), 1 + j * n)
    Z = Vector{Float64}(undef, 1 + j * n)
    # Step 1: Initialize
    s = deepcopy(sₒ)
    s_best = deepcopy(sₒ)
    x = h(s)
    z = f(s)
    X[1] = x
    Z[1] = z
    z_best = z
    t = ω̅ * f(sₒ) / log(1 / τ̅)
    Cᵣ = zeros(Int, Rᵣ, Rᵪ)
    Cᵢ = zeros(Int, Iᵣ, Iᵪ)
    Pᵣ = zeros(Float64, Rᵣ, Rᵪ)
    Pᵢ = zeros(Float64, Iᵣ, Iᵪ)
    Wᵣ = ones(Float64, Rᵣ, Rᵪ)
    Wᵢ = ones(Float64, Iᵣ, Iᵪ)
    Sᵣ = zeros(Float64, Rᵣ, Rᵪ)
    Sᵢ = zeros(Float64, Iᵣ, Iᵪ)
    # Step 2: Loop over segments
    for segment ∈ 1:j
        # Step 2.1: Reset count and score for every removal and insertion operator
        Cᵣ .= 0; Sᵣ .= 0
        Cᵢ .= 0; Sᵢ .= 0
        # Step 2.2: Update selection probability for every removal and insertion operator
        for row ∈ 1:Rᵣ
            row_sum = sum(Wᵣ[row, :])
            for col ∈ 1:Rᵪ
                Pᵣ[row, col] = Wᵣ[row, col] / row_sum
            end
        end
        for row ∈ 1:Iᵣ
            row_sum = sum(Wᵢ[row, :])
            for col ∈ 1:Iᵪ
                Pᵢ[row, col] = Wᵢ[row, col] / row_sum
            end
        end

        # Step 2.3: Loop over iterations within the segment
        for iteration ∈ 1:n
            # Step 2.3.1: Randomly select a removal and an insertion operator based on operator selection probabilities, and consequently update count for the selected operators
            r = [sample(rng, 1:Rᵪ, Weights(Pᵣ[row, :])) for row ∈ 1:Rᵣ]
            i = [sample(rng, 1:Iᵪ, Weights(Pᵢ[row, :])) for row ∈ 1:Iᵣ]
            for row ∈ 1:Rᵣ
                Cᵣ[row, r[row]] += 1
            end
            for row ∈ 1:Iᵣ
                Cᵢ[row, i[row]] += 1
            end
            # Step 2.3.2: Using the selected removal and insertion operators destroy and repair the current solution to develop a new solution
            q = removalcount(rng, s, μ̲, e̲, μ̅, e̅)
            s′ = deepcopy(s)
            for row ∈ 1:Rᵣ
                s′ = Ψᵣ[row, r[row]](rng, q, s′)
            end
            for row ∈ 1:Iᵣ
                s′ = Ψᵢ[row, i[row]](rng, s′)
            end
            x′ = h(s′)
            z′ = f(s′)
            # Step 2.3.3: If this new solution is better than the best solution, then set the best solution and the current solution to the new solution, and accordingly update scores of the selected removal and insertion operators by σ₁
            if z′ < z_best
                s_best = s′
                s = s′
                x = x′
                z = z′
                z_best = z′
                for row ∈ 1:Rᵣ
                    Sᵣ[row, r[row]] += σ₁
                end
                for row ∈ 1:Iᵣ
                    Sᵢ[row, i[row]] += σ₁
                end
            # Step 2.3.4: Else if this new solution is only better than the current solution, then set the current solution to the new solution and accordingly update scores of the selected removal and insertion operators by σ₂
            elseif z′ < z
                s = s′
                x = x′
                z = z′
                if x′ ∉ X
                    for row ∈ 1:Rᵣ
                        Sᵣ[row, r[row]] += σ₂
                    end
                    for row ∈ 1:Iᵣ
                        Sᵢ[row, i[row]] += σ₂
                    end
                end
            # Step 2.3.5: Else accept the new solution with simulated annealing acceptance criterion. Further, if the new solution is also newly found then update operator scores by σ₃
            else 
                if rand(rng) < exp((z - z′) / t)
                    s = s′
                    x = x′
                    z = z′
                    if x′ ∉ X
                        for row ∈ 1:Rᵣ
                            Sᵣ[row, r[row]] += σ₃
                        end
                        for row ∈ 1:Iᵣ
                            Sᵢ[row, i[row]] += σ₃
                        end
                    end
                end
            end
            # Step 2.3.6: Update temperature according to cooling schedule
            t = max(ω̲  * f(s_best) / log(1 / τ̲), t * θ) 
            # Step 2.3.7: Store best solution from each iteration
            X[1 + (segment - 1) * n + iteration] = x
            Z[1 + (segment - 1) * n + iteration] = z
        end
        # Step 2.4: Update weights for every removal and insertion operator (module)
        for row ∈ 1:Rᵣ
            for col ∈ 1:Rᵪ
                if !iszero(Cᵣ[row, col])
                    Wᵣ[row, col] = (1 - ρ) * Wᵣ[row, col] + ρ * Sᵣ[row, col] / Cᵣ[row, col]
                end
            end
        end
        for row ∈ 1:Iᵣ
            for col ∈ 1:Iᵪ
                if !iszero(Cᵢ[row, col])
                    Wᵢ[row, col] = (1 - ρ) * Wᵢ[row, col] + ρ * Sᵢ[row, col] / Cᵢ[row, col]
                end
            end
        end
        # Step 2.5: Reset current solution
        s = sₒ
        # Step 2.6: Perform local search
        s′ = deepcopy(s)
        for l ∈ L
            s′ = Ψₗ[l](rng, m, s′)
        end
        x′ = h(s′)
        z′ = f(s′)
        if z′ < z_best
            s_best = s′
            s = s′
            x = x′
            z = z′
            z_best = z′
        elseif z′ < z
            s = s′
            x = x′
            z = z′
        end
    end
    # Step 3: Display the convergence plot and return the best solution
    return s_best
end
