"""
    ALNS([rng::AbstractRNG], χ::ALNSparameters, sₒ::Solution; mute=false))

Adaptive Large Neighborhood Search (ALNS)

Given ALNS optimization parameters `χ` and an initial solution `sₒ`, 
ALNS adaptively searches large neighborhoods in the solution domain and
returns the best found solution. Additionally, displays a convergence plot.

Takes `mute` as argument. If `true` mutes progressbar and pltcnv output.

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
function ALNS(rng::AbstractRNG, χ::ALNSparameters, sₒ::Solution; mute=false)
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
    # Step 1: Initialize
    s = sₒ
    s⃰ = sₒ
    t = ω̅  * Δ/log(1/τ̅)
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
            r = rand(rng, R, Weights(Pᵣ))
            i = rand(rng, I, Weights(Pᵢ))
            Cᵣ[r] += 1
            Cᵢ[i] += 1
            # Step 2.3.2: Using the selected removal and insertion operators destroy and repair the current solution to develop a new solution
            
            # Step 2.3.3: If this new solution is better than the best solution, then set the best solution and the current solution to the new solution, and accordingly update scores of the selected removal and insertion operators by σ₁
            if f(s') < f(s⃰)
                s⃰ = s′
                s = s′
                Sᵣ[r] += σ₁
                Sᵢ[i] += σ₁
            # Step 2.3.4: Else if this new solution is only better than the current solution, then set the current solution to the new solution and accordingly update scores of the selected removal and insertion operators by σ₂
            elseif f(s′) < f(s)
                s = s′
                Sᵣ[r] += σ₂
                Sᵢ[i] += σ₂
            # Step 2.3.5: Else accept the new solution with simulated annealing acceptance criterion. Further, if the new solution is also newly found then update operator scores by σ₃
            else 
                if rand(rng) < exp((f(s) - f(s′)) / t)
                    s = s′
                    Sᵣ[r] += σ₃
                    Sᵢ[i] += σ₃
                end
            end

        end
        # Step 2.4: Update weights for every removal and insertion operator
        for r ∈ R if !iszero(Cᵣ[r]) Wᵣ[r] = (1 - ρ) * Wᵣ[r] + ρ * Sᵣ[r] / Cᵣ[r] end end
        for i ∈ I if !iszero(Cᵢ[i]) Wᵢ[i] = (1 - ρ) * Wᵢ[i] + ρ * Sᵢ[i] / Cᵢ[i] end end
        # Step 2.5: Reset current solution
        s = sₒ
        # Step 2.6: Perform local search
        s' = s
        for l ∈ L """do local search""" end
        if f(s′) < f(s⃰)
            s⃰ = s′
            s = s′
        elseif f(s′) < f(s)
            s = s′
        end
    end
    # Step 3: Display the convergence plot and return the best solution
    return s⃰
end