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
    # TODO: Implement ALNS here
    # Step 0: Pre-initialize
    # Step 1: Initialize
    # Step 2: Loop over segments
    # Step 2.1: Reset count and score for every removal and insertion operator
    # Step 2.2: Update selection probability for every removal and insertion operator
    # Step 2.3: Loop over iterations within the segment
    # Step 2.3.1: Randomly select a removal and an insertion operator based on operator selection probabilities, and consequently update count for the selected operators
    # Step 2.3.2: Using the selected removal and insertion operators destroy and repair the current solution to develop a new solution
    # Step 2.3.3: If this new solution is better than the best solution, then set the best solution and the current solution to the new solution, and accordingly update scores of the selected removal and insertion operators by σ₁
    # Step 2.3.4: Else if this new solution is only better than the current solution, then set the current solution to the new solution and accordingly update scores of the selected removal and insertion operators by σ₂
    # Step 2.3.5: Else accept the new solution with simulated annealing acceptance criterion. Further, if the new solution is also newly found then update operator scores by σ₃
    # Step 2.4: Update weights for every removal and insertion operator
    # Step 2.5: Reset current solution
    # Step 2.6: Perform local search
    # Step 3: Display the convergence plot and return the best solution
    return s⃰
end