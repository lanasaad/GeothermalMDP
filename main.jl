using Distributions
using POMDPTools
using Plots
using Random
using MCTS

include("src/pomdp.jl")
include("src/policy.jl")
include("src/utils.jl")

rng = MersenneTwister(1)
config = GeothermalConfig(Δrow=100)
mdp = GeothermalMDP(df=config.df)
mdp.df
# Random Policy
policy_random = BuildRandom(mdp)
hr_random = HistoryRecorder(max_steps = 30)
rhist_random = simulate(hr_random, mdp, policy_random)

# Track rewards during the simulation
rewards = Float64[]

for step in rhist_random
    reward = POMDPs.reward(mdp, step.s, step.a)
    push!(rewards, reward)
end

# Analyze NPV, mean value, and standard deviation
npv_random = sum(rewards)
mean_random = mean(rewards)
stddev_random = std(rewards)

println("Random Policy - NPV: $npv_random, Mean: $mean_random, StdDev: $stddev_random")

# Save GIF
filename = "figs/geothermal_random.gif"
save_gif(rhist_random, filename, 2024, 4)
println("Done, GIF saved to $filename.")


# Define the MCTS solver
solver_mcts = MCTSSolver(n_iterations=5, depth=3, rng=rng)

# Solve for the policy using MCTS
policy_mcts = solve(solver_mcts, mdp)

# Simulate the policy
hr_mcts = HistoryRecorder(max_steps = 30)
rhist_mcts = simulate(hr_mcts, mdp, policy_mcts)

# Save the simulation as a GIF
filename_mcts = "figs/geothermal_mcts.gif"
save_gif(rhist_mcts, filename_mcts, 2024, 4)
println("Done, MCTS Policy GIF saved to $filename_mcts.")