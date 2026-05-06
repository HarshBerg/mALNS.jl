using mALNS
using CSV
using Revise
using Random
using DataFrames

let
    # Define instances
    instances = ["X-n411-k19", "X-n420-k130", "X-n439-k37", "X-n469-k138", "X-n480-k70", "X-n502-k39", "X-n524-k153", "X-n548-k50", "X-n561-k42", "X-n586-k159",
                 "X-n613-k62", "X-n627-k43", "X-n641-k35", "X-n655-k131", "X-n670-k130", "X-n701-k44", "X-n716-k35", "X-n733-k159", "X-n749-k98", "X-n766-k71",
                 "X-n801-k40", "X-n819-k171", "X-n856-k95", "X-n876-k59", "X-n895-k37", "X-n916-k207", "X-n936-k151", "X-n957-k87", "X-n979-k58", "X-n1001-k43"]
    # Define random number generator seeds
    seeds = [1234, 1729, 2310, 3103, 9999, 4200, 5544, 7788, 8080, 6000]
    # Dataframes to store solution quality and run time
    dfᶠ = DataFrame([instances, zeros(length(instances)), [zeros(length(instances)) for _ ∈ seeds]...], ["instance", "initial", ["$seed" for seed ∈ seeds]...])
    dfᵗ = DataFrame([instances, zeros(length(instances)), [zeros(length(instances)) for _ ∈ seeds]...], ["instance", "initial", ["$seed" for seed ∈ seeds]...])
    for (i,instance) ∈ enumerate(instances)
        # Visualize instance
        display(visualize(instance))
        # Define inital solution method and build the initial solution
        t₁ = @elapsed s₁ = initialize(instance);
        # Visualize initial solution
        display(visualize(s₁))
        # Fetch solution characteristics
        println("\nINSTANCE: $instance")
        println("Initial Solution:")
        println("   Run Time: $(round(t₁, digits=3)) seconds")
        println("   Feasiblity: $(isfeasible(s₁)) | $(round(s₁.p, digits=3))")
        println("   Objective Function: $(round(f(s₁), digits=3))")
        # Store results
        dfᶠ[i,2] = f(s₁)
        dfᵗ[i,2] = t₁
        # Save results
        CSV.write("objective_function.csv", dfᶠ)
        CSV.write("run_time.csv", dfᵗ)
        for (j,seed) ∈ enumerate(seeds)
            println("\nOptimal Solution | seed: $seed")
            rng = MersenneTwister(seed);
            # Define ALNS parameters
            x = max(100, lastindex(s₁.G.N))
            χ = ALNSparameters(
                j   =   100                       ,
                k   =   10                      ,
                n   =   x                       ,
                m   =   100x                    ,
                Ψᵣ  =   [
                            randomnode!         randomarc!          randomsegment! ;
                            relatednode!        relatedarc!         relatedsegment!;
                            worstnode!          worstarc!           worstsegment!
                        ]                       ,
                Ψᵢ  =   [
                            bestprecise!        bestperturb!        ;
                            greedyprecise!      greedyperturb!      ;
                            regret2precise!     regret2perturb!     ;
                            regret3precise!     regret3perturb!     ;
                            regret5precise!     regret5perturb!     
                        ]                       ,
                Ψₗ  =   [
                            intramove!          ,
                            intermove!          ,
                            intraswap!          ,
                            interswap!          ,
                            intraopt!           ,
                            interopt!           ,
                        ]                       ,
                σ₁  =   15.0                    ,
                σ₂  =   10.0                    ,
                σ₃  =   3.0                    ,
                μ̲   =   0.05                    ,
                e̲   =   5                       ,
                μ̅   =   0.3                     ,
                e̅   =   30                      ,
                ω̅   =   1/x                     ,
                τ̅   =   0.1                     ,
                ω̲   =   0.1/x                   ,
                τ̲   =   0.01                    ,
                θ   =   exp(log(0.1 / 1 * log(1 / 0.1) / log(1 / 0.01)) / (x * x))                  ,
                ρ   =   0.1
            );
            # Run ALNS and fetch best solution
            t₂ = @elapsed s₂ = conALNS(rng, χ, s₁);
            # Visualize best solution
            display(visualize(s₂))
            # Fetch solution characteristics
            println("   Run Time: $(round(t₂, digits=3)) seconds")
            println("   Feasiblity: $(isfeasible(s₂)) | $(round(s₂.p, digits=3))")
            println("   Objective Function: $(round(f(s₂), digits=3))")
            # Store results
            dfᶠ[i,j+2] = f(s₂)
            dfᵗ[i,j+2] = t₂
            # Save results
            # sol(s₂, "solutions/$instance-seed$seed")
            CSV.write("objective_function.csv", dfᶠ)
            CSV.write("run_time.csv", dfᵗ)
        end
        println(dfᶠ)
        println(dfᵗ)
    end
    return
end 