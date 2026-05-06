"""
    conALNS([rng::AbstractRNG], χ::ALNSparameters, sₒ::Solution; mute=false))

conventional Adaptive Large Neighborhood Search (conALNS)

Given ALNS optimization parameters `χ` and an initial solution `sₒ`, 
ALNS adaptively searches large neighborhoods in the solution domain and
returns the best found solution. Additionally, displays a convergence plot.

Takes `mute` as argument. If `true` mutes progressbar and pltcnv output.

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
function conALNS(rng::AbstractRNG, χ::ALNSparameters, sₒ::Solution; mute=false)
    # Step 0: Pre-initialize
    G = sₒ.G
    j, k = χ.j, χ.k
    n, m = χ.n, χ.m
    Ψᵣ, Ψᵢ, Ψₗ = χ.Ψᵣ, χ.Ψᵢ, χ.Ψₗ
    σ₁, σ₂, σ₃ = χ.σ₁, χ.σ₂, χ.σ₃
    μ̲, e̲ = χ.μ̲, χ.e̲
    μ̅, e̅ = χ.μ̅, χ.e̅
    ω̅, τ̅ = χ.ω̅, χ.τ̅
    ω̲, τ̲ = χ.ω̲, χ.τ̲
    θ, ρ = χ.θ, χ.ρ
    e = lastindex(G.N)
    R = eachindex(Ψᵣ)
    I = eachindex(Ψᵢ)
    L = eachindex(Ψₗ)
    X = Vector{UInt}(undef, 1 + j * (n + 1))
    Z = Vector{Float64}(undef, 1 + j * (n + 1))
    G = sₒ.G
    # Step 1: Initialize
    s = deepcopy(sₒ)
    x = h(s)
    z = f(s)
    X[1] = x
    Z[1] = z
    s⃰ = s
    z⃰ = z
    t = ω̅ * z⃰/log(1/τ̅)
    Cᵣ, Pᵣ, Πᵣ, Wᵣ = zeros(Int, R), zeros(R), zeros(R), ones(R)
    Cᵢ, Pᵢ, Πᵢ, Wᵢ = zeros(Int, I), zeros(I), zeros(I), ones(I)
    φ = 0
    # Step 2: Loop over segments.
    if !mute p = Progress(n * j, desc="Computing...", color=:blue, showspeed=true) end
    for u ∈ 1:j
        # Step 2.1: Reset count and score for every removal and insertion operator
        for r ∈ R Cᵣ[r], Πᵣ[r] = 0, 0. end
        for i ∈ I Cᵢ[i], Πᵢ[i] = 0, 0. end
        # Step 2.2: Update selection probability for every removal and insertion operator
        for r ∈ R Pᵣ[r] = Wᵣ[r] / sum(Wᵣ) end
        for i ∈ I Pᵢ[i] = Wᵢ[i] / sum(Wᵢ) end
        # Step 2.3: Loop over iterations within the segment
        for v ∈ 1:n
            # Step 2.3.1: Randomly select a removal and an insertion operator based on operator selection probabilities, and consequently update count for the selected operators.
            r = sample(rng, R, Weights(Pᵣ))
            i = sample(rng, I, Weights(Pᵢ))
            Cᵣ[r] += 1
            Cᵢ[i] += 1
            # Step 2.3.2: Using the selected removal and insertion operators destroy and repair the current solution to develop a new solution.
            η = rand(rng)
            q = Int(floor(((1 - η) * min(e̲, μ̲ * e) + η * min(e̅, μ̅ * e))))
            s′= deepcopy(s)
            remove!(rng, q, s′, Ψᵣ[r])
            insert!(rng, s′, Ψᵢ[i])
            x′ = h(s′)
            z′ = f(s′)
            # Step 2.3.3: If this new solution is better than the best solution, then set the best solution and the current solution to the new solution, and accordingly update scores of the selected removal and insertion operators by σ₁.
            if z′ < z⃰
                s = s′
                s⃰ = s′
                z = z′
                z⃰ = z′
                Πᵣ[r] += σ₁
                Πᵢ[i] += σ₁
                φ = 1
            # Step 2.3.4: Else if this new solution is only better than the current solution, then set the current solution to the new solution and accordingly update scores of the selected removal and insertion operators by σ₂.
            elseif z′ < z
                s = s′
                z = z′
                if x′ ∉ X
                    Πᵣ[r] += σ₂
                    Πᵢ[i] += σ₂
                end
                φ = φ
            # Step 2.3.5: Else accept the new solution with simulated annealing acceptance criterion. Further, if the new solution is also newly found then update operator scores by σ₃.
            else
                η = rand(rng)
                if η < exp(-(z′ - z)/t)
                    s = s′
                    z = z′
                    if x′ ∉ X
                        Πᵣ[r] += σ₃
                        Πᵢ[i] += σ₃
                    end
                    φ = φ
                end
            end
            x = x′
            X[1 + (u - 1) * (n + 1) + v] = x
            Z[1 + (u - 1) * (n + 1) + v] = z
            t = max(t * θ, ω̲ * z⃰/log(1/τ̲))
            if !mute next!(p) end
        end
        # Step 2.4: Update weights for every removal and insertion operator.
        for r ∈ R if !iszero(Cᵣ[r]) Wᵣ[r] = ρ * Πᵣ[r] / Cᵣ[r] + (1 - ρ) * Wᵣ[r] end end
        for i ∈ I if !iszero(Cᵢ[i]) Wᵢ[i] = ρ * Πᵢ[i] / Cᵢ[i] + (1 - ρ) * Wᵢ[i] end end
        # Step 2.5: Reset current solution.
        if iszero(u % k)
            s = iszero(φ) ? deepcopy(s⃰) : s
            z = iszero(φ) ? z⃰ : z
            φ = 0
        end
        # Step 2.6: Perform local search.
        s′ = deepcopy(s)
        for l ∈ L localsearch!(rng, m, s′, Ψₗ[l]) end
        x′ = h(s′)
        z′ = f(s′)
        if z′ < z⃰
            s = s′
            s⃰ = s′
            z = z′
            z⃰ = z′
        elseif z′ < z
            s = s′
            z = z′
        end
        x = x′
        X[1 + u * (n + 1)] = x
        Z[1 + u * (n + 1)] = z
    end
    # Step 3: Display the convergence plot and return the best solution
    if !mute display(pltcnv(X,Z)) end
    return s⃰
end
conALNS(χ::ALNSparameters, s::Solution; mute=false) = conALNS(Random.GLOBAL_RNG, χ, s; mute=mute)

"""
    modALNS([rng::AbstractRNG], χ::ALNSparameters, sₒ::Solution; mute=false))

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
    G = sₒ.G
    j, k = χ.j, χ.k
    n, m = χ.n, χ.m
    Ψᵣ, Ψᵢ, Ψₗ = χ.Ψᵣ, χ.Ψᵢ, χ.Ψₗ
    σ₁, σ₂, σ₃ = χ.σ₁, χ.σ₂, χ.σ₃
    μ̲, e̲ = χ.μ̲, χ.e̲
    μ̅, e̅ = χ.μ̅, χ.e̅
    ω̅, τ̅ = χ.ω̅, χ.τ̅
    ω̲, τ̲ = χ.ω̲, χ.τ̲
    θ, ρ = χ.θ, χ.ρ
    e = lastindex(G.N)
    rₚ, rₛ = size(Ψᵣ)
    iₚ, iₛ = size(Ψᵢ)
    L = eachindex(Ψₗ)
    X = Vector{UInt}(undef, 1 + j * (n + 1))
    Z = Vector{Float64}(undef, 1 + j * (n + 1))
    G = sₒ.G
    # Step 1: Initialize
    s = deepcopy(sₒ)
    x = h(s)
    z = f(s)
    X[1] = x
    Z[1] = z
    s⃰ = s
    z⃰ = z
    t = ω̅ * z⃰/log(1/τ̅)
    Cᵣₚ, Pᵣₚ, Πᵣₚ, Wᵣₚ = zeros(Int, rₚ), zeros(rₚ), zeros(rₚ), ones(rₚ)
    Cᵣₛ, Pᵣₛ, Πᵣₛ, Wᵣₛ = zeros(Int, rₛ), zeros(rₛ), zeros(rₛ), ones(rₛ)
    Cᵢₚ, Pᵢₚ, Πᵢₚ, Wᵢₚ = zeros(Int, iₚ), zeros(iₚ), zeros(iₚ), ones(iₚ)
    Cᵢₛ, Pᵢₛ, Πᵢₛ, Wᵢₛ = zeros(Int, iₛ), zeros(iₛ), zeros(iₛ), ones(iₛ)
    φ = 0
    # Step 2: Loop over segments.
    if !mute p = Progress(n * j, desc="Computing...", color=:blue, showspeed=true) end
    for u ∈ 1:j
        # Step 2.1: Reset count and score for every removal and insertion operator
        Cᵣₚ .= 0; Πᵣₚ .= 0
        Cᵢₚ .= 0; Πᵢₚ .= 0
        # Step 2.2: Update selection probability for every removal and insertion operator
        Pᵣₚ = Wᵣₚ / sum(Wᵣₚ)
        Pᵣₛ = Wᵣₛ / sum(Wᵣₛ)
        Pᵢₚ = Wᵢₚ / sum(Wᵢₚ)
        Pᵢₛ = Wᵢₛ / sum(Wᵢₛ)
        # Step 2.3: Loop over iterations within the segment
        for v ∈ 1:n
            # Step 2.3.1: Randomly select a removal and an insertion operator based on operator selection probabilities, and consequently update count for the selected operators
            pᵣ, sᵣ = sample(rng, 1:rₚ, Weights(Pᵣₚ)), sample(rng, 1:rₛ, Weights(Pᵣₛ))
            pᵢ, sᵢ = sample(rng, 1:iₚ, Weights(Pᵢₚ)), sample(rng, 1:iₛ, Weights(Pᵢₛ))
            # Step 2.3.2: Using the selected removal and insertion operators destroy and repair the current solution to develop a new solution
            η = rand(rng)
            q = Int(floor((1 - η) * min(e̲, μ̲  * e) + η * min(e̅, μ̅  * e)))
            s′ = deepcopy(s)
            remove!(rng, q, s′, Ψᵣ[pᵣ,sᵣ])
            insert!(rng, s′, Ψᵢ[pᵢ,sᵢ])
            x′ = h(s′)
            z′ = f(s′)
            # Step 2.3.3: If this new solution is better than the best solution, then set the best solution and the current solution to the new solution, and accordingly update scores of the selected removal and insertion operators by σ₁.
            if z′ < z⃰
                s = s′
                s⃰ = s′
                z = z′
                z⃰ = z′
                Πᵣₚ[pᵣ] += σ₁
                Πᵣₛ[sᵣ] += σ₁
                Πᵢₚ[pᵢ] += σ₁
                Πᵢₛ[sᵢ] += σ₁
                Cᵣₚ[pᵣ] += 1
                Cᵣₛ[sᵣ] += 1
                Cᵢₚ[pᵢ] += 1
                Cᵢₛ[sᵢ] += 1
                φ = 1
            # Step 2.3.4: Else if this new solution is only better than the current solution, then set the current solution to the new solution and accordingly update scores of the selected removal and insertion operators by σ₂.
            elseif z′ < z
                s = s′
                z = z′
                if x′ ∉ X
                    Πᵣₚ[pᵣ] += σ₂
                    Πᵣₛ[sᵣ] += σ₂
                    Πᵢₚ[pᵢ] += σ₂
                    Πᵢₛ[sᵢ] += σ₂
                    Cᵣₚ[pᵣ] += 1
                    Cᵣₛ[sᵣ] += 1
                    Cᵢₚ[pᵢ] += 1
                    Cᵢₛ[sᵢ] += 1
                end
                φ = φ
            # Step 2.3.5: Else accept the new solution with simulated annealing acceptance criterion. Further, if the new solution is also newly found then update operator scores by σ₃.
            else
                η = rand(rng)
                if η < exp(-(z′ - z)/t)
                    s = s′
                    z = z′
                    if x′ ∉ X
                        Πᵣₚ[pᵣ] += σ₃
                        Πᵣₛ[sᵣ] += σ₃
                        Πᵢₚ[pᵢ] += σ₃
                        Πᵢₛ[sᵢ] += σ₃
                        Cᵣₚ[pᵣ] += 1
                        Cᵣₛ[sᵣ] += 1
                        Cᵢₚ[pᵢ] += 1
                        Cᵢₛ[sᵢ] += 1
                    end
                    φ = φ
                end
            end
            x = x′
            X[1 + (u - 1) * (n + 1) + v] = x
            Z[1 + (u - 1) * (n + 1) + v] = z
            t = max(t * θ, ω̲ * z⃰/log(1/τ̲))
            if !mute next!(p) end
        end
        # Step 2.4: Update weights for every removal and insertion operator (module)
        for pᵣ ∈ 1:rₚ Wᵣₚ[pᵣ] = iszero(Cᵣₚ[pᵣ]) ? Wᵣₚ[pᵣ] : (1 - ρ) * Wᵣₚ[pᵣ] + ρ * Πᵣₚ[pᵣ] / Cᵣₚ[pᵣ] end
        for sᵣ ∈ 1:rₛ Wᵣₛ[sᵣ] = iszero(Cᵣₛ[sᵣ]) ? Wᵣₛ[sᵣ] : (1 - ρ) * Wᵣₛ[sᵣ] + ρ * Πᵣₛ[sᵣ] / Cᵣₛ[sᵣ] end
        for pᵢ ∈ 1:iₚ Wᵢₚ[pᵢ] = iszero(Cᵢₚ[pᵢ]) ? Wᵢₚ[pᵢ] : (1 - ρ) * Wᵢₚ[pᵢ] + ρ * Πᵢₚ[pᵢ] / Cᵢₚ[pᵢ] end
        for sᵢ ∈ 1:iₛ Wᵢₛ[sᵢ] = iszero(Cᵢₛ[sᵢ]) ? Wᵢₛ[sᵢ] : (1 - ρ) * Wᵢₛ[sᵢ] + ρ * Πᵢₛ[sᵢ] / Cᵢₛ[sᵢ] end
        # Step 2.5: Reset current solution.
        if iszero(u % k)
            s = iszero(φ) ? deepcopy(s⃰) : s
            z = iszero(φ) ? z⃰ : z
            φ = 0
        end
        # Step 2.6: Perform local search.
        s′ = deepcopy(s)
        for l ∈ L localsearch!(rng, m, s′, Ψₗ[l]) end
        x′ = h(s′)
        z′ = f(s′)
        if z′ < z⃰
            s = s′
            s⃰ = s′
            z = z′
            z⃰ = z′
        elseif z′ < z
            s = s′
            z = z′
        end
        x = x′
        X[1 + u * (n + 1)] = x
        Z[1 + u * (n + 1)] = z
    end
    # print probalities of primary and secondary modules for removal and insertion operators
    #println("Removal Operator Probabilities (Primary): ", Pᵣₚ)
    #println("Removal Operator Probabilities (Secondary): ", Pᵣₛ)
    #println("Insertion Operator Probabilities (Primary): ", Pᵢₚ)
    #println("Insertion Operator Probabilities (Secondary): ", Pᵢₛ)
    # Step 3: Display the convergence plot and return the best solution
    if !mute display(pltcnv(X,Z)) end
    return s⃰
end
modALNS(χ::ALNSparameters, s::Solution; mute=false) = modALNS(Random.GLOBAL_RNG, χ, s; mute=mute)
