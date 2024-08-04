using Distributions
using POMDPTools
using Plots
using Random

include("src/pomdp.jl")
include("src/policy.jl")
include("src/utils.jl")

rng = MersenneTwister(1)
mdp = GeothermalMDP()

# random policy
policy = BuildRandom(mdp)
hr = HistoryRecorder(max_steps = 30)
@time rhist = simulate(hr, mdp, policy);

# save gif
filename = "figs/geothermal.gif"
year = 2024
fps = 4
save_gif(rhist, filename, year, fps)

println("Done, gif saved to $filename.")