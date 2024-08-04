using POMDPs
using Parameters
using Random
using CSV
using DataFrames
using Revise
using IterTools
using DifferentialEquations
using DataFrames
using Distances
using Distributions
using POMDPTools
using Plots


@with_kw mutable struct GeothermalConfig
	inputfile::String = "data/reservoir.csv"
    Δrow = 10
    data = CSV.read(inputfile, DataFrame)
    df = data[1:Δrow:end, :]
end


@with_kw mutable struct Coord
	x::Float64
	y::Float64
end

@with_kw mutable struct State
	t::Vector{Float64}
	i::Vector{Coord}
    p::Vector{Coord}
end


@with_kw mutable struct DrillAction
	i::Coord
	p::Coord
end


@with_kw mutable struct ShutAction
	i::Coord
	p::Coord
end


@with_kw mutable struct NothingAction
	i::Coord = Coord(-1, -1)
	p::Coord = Coord(-1, -1)
end


const Action = Union{DrillAction, ShutAction, NothingAction}


#specify injector and producer locations here to populate state space
@with_kw mutable struct GeothermalMDP <: MDP{State, Action}
    df::DataFrame = GeothermalConfig().df   #dataframe
    xs::Vector{Float64} = df.x              #x coordinates
    ys::Vector{Float64} = df.y              #y coordinates
    t0::Vector{Float64} = df.temperature    #initial temperature
    w0::Vector{Coord} = []                  #initial well locations
    t_values::Vector{Float64} = 30:0.01:140 #temperature values
    num_max_well_pairs::Int = 99            #maximum number of well pairs
    min_well_distance::Float64 = 5.0        #minimum distance between wells
	capex::Float64 = 10000                  #capital expenditure for each well pair
	opex::Float64 = 1000                    #operational expenditure for each well
    μ_Δtemp::Float64 = 1.0                  #mean base temperature change if nearby wells are present
    σ_Δtemp::Float64 = 1.0                  #standard deviation of temperature change if nearby wells are present
    d_Δtemp::Float64 = 10.0                 #distances affected by nearby wells
    dissipation::Float64 = 0.95             #temperature dissipation factor
	revenue::Float64 = 5000                 #revenue per well pair, if temperature ≥ opt_temp
    opt_temp::Float64 = 80.0                #optimal temperature for revenue
    Δrev::Float64 = 100                     #revenue change per degree below opt_temp
	γ::Float64 = 0.95                       #discount factor
    coords = [Coord(x, y) for (x,y) in zip(xs, ys)]
    well_pairs = [(i, j) for (i, j) in IterTools.product(coords, coords) if i != j]
    rng = MersenneTwister(2024)             #random number generator
end



function POMDPs.states(mdp::GeothermalMDP)
	return nothing #TODO: add if needed by online solvers
end


function POMDPs.actions(mdp::GeothermalMDP)
    drill_actions = [DrillAction(i=coord_1, p=coord_2) for (coord_1, coord_2) in mdp.well_pairs]
    shut_actions = [ShutAction(i=coord_1, p=coord_2) for (coord_1, coord_2) in mdp.well_pairs]
    nothing_actions = [NothingAction()]
    actions = vcat(drill_actions, shut_actions, nothing_actions)
	return actions
end


function POMDPs.actions(mdp::GeothermalMDP, s::State)
    well_pairs = deepcopy(mdp.well_pairs)

    #remove well pairs where either injector or producer is already present
    for (i, j) in well_pairs
        if Coord(i.x, i.y) in s.i || Coord(j.x, j.y) in s.p
            deleteat!(well_pairs, findall(x -> x == (i, j), well_pairs))
        end
    end

    drill_actions = [DrillAction(i=coord_1, p=coord_2) for (coord_1, coord_2) in well_pairs]

    nothing_actions = [NothingAction()]

    #shut actions only at locations where there are wells
    if length(s.i) > 0 && length(s.p) > 0
        shut_actions = [[ShutAction(i=coord_1, p=coord_2) for (coord_1, coord_2) in IterTools.product(s.i, s.p)]...]
        all_actions = vcat(drill_actions, shut_actions, nothing_actions)
    else
        shut_actions = []
        all_actions = vcat(drill_actions, nothing_actions)
    end

    return all_actions
end

function p(mdp::GeothermalMDP, well::Coord, coord::Coord, λ::Float64=5e2)
    dist_ = euclidean([coord.x, coord.y], [well.x, well.y])
    return (coord ∈ neighbors(mdp, well)) * (1 - exp(-λ/dist_))
end


coord2index(mdp::GeothermalMDP, coord::Coord) = findfirst(x -> isapprox(x.x, coord.x) && isapprox(x.y, coord.y), mdp.coords)


function neighbors(mdp::GeothermalMDP, well::Coord, d::Float64=500.0)
    return [coord for coord in mdp.coords if euclidean([coord.x, coord.y], [well.x, well.y]) <= d]
end


function update_neighbor_temp(t0::Vector{Float64}, mdp::GeothermalMDP, w::Coord, Δt::Float64, rng::AbstractRNG)
    t1 = deepcopy(t0)
    for n in neighbors(mdp, w)
        if rand(rng) < p(mdp, w, n, 400.)
            t1[coord2index(mdp, n)] += Δt * p(mdp, w, n)
        end
    end
    return t1
end

function step_temp(mdp::GeothermalMDP, injs::Vector{Coord}, prods::Vector{Coord}, temps::Vector{Float64}, rng::AbstractRNG)
    #TODO: move to MDP struct
    Δtemp = 3.0
    next_temp = deepcopy(temps)
    
    #process the effect of injection wells
    for inj in injs
        next_temp = update_neighbor_temp(next_temp, mdp, inj, -Δtemp, rng)
    end
    #TODO: add the effect of production wells
    #TODO: add natural dissipation
    return next_temp
end

POMDPs.isterminal(mdp::GeothermalMDP, s::State) = length(s.i) >= mdp.num_max_well_pairs || length(s.p) >= mdp.num_max_well_pairs

function POMDPs.transition(mdp::GeothermalMDP, s::State, a::Action)

    rng = mdp.rng
	if isterminal(mdp, s)
        return Deterministic(s)
    end

	if typeof(a) == DrillAction
		injectors = vcat(s.i, [a.i])
        producers = vcat(s.p, [a.p])
        next_temp = step_temp(mdp, injectors, producers, s.t, rng)
        # println("DrillAction: ", a.i, " ", a.p)
        return Deterministic(State(t = next_temp, i = injectors, p = producers)) #step_temp is already random

    elseif typeof(a) == ShutAction
        injectors = [i for i in s.i if i != a.i]
        producers = [p for p in s.p if p != a.p]
        next_temp = step_temp(mdp, injectors, producers, s.t, rng)
        # println("ShutAction: ", a.i, " ", a.p)
        return Deterministic(State(t = next_temp, i = injectors, p = producers))

    else #NothingAction
        next_temp = step_temp(mdp, s.i, s.p, s.t, rng)
        # println("NothingAction")
        return Deterministic(State(t = next_temp, i = s.i, p = s.p))
    end

end


function compute_revenue(mdp::GeothermalMDP, temp::Vector{Float64}, producers::Vector{Coord})

    revenue = 0.0
    for p in producers
        idx = coord2index(mdp, p)
        if isnothing(idx)
            println("Error: Index not found for ", p)

        elseif temp[idx] >= mdp.opt_temp
            revenue += mdp.revenue
        else
            mdp.Δrev = (mdp.opt_temp - temp[idx]) / mdp.σ_Δtemp
        end          
    end
    return revenue
end


function POMDPs.reward(mdp::GeothermalMDP, s::State, a::Action)
    if isterminal(mdp, s)
        return 0.0
    end

    if typeof(a) == DrillAction
        producers = vcat(s.p, [a.p])
    elseif typeof(a) == ShutAction
        producers = [p for p in s.p if p != a.p]
    else
        producers = s.p
    end

    revenue = compute_revenue(mdp, s.t, producers)
    capex = typeof(a) == DrillAction ? mdp.capex : 0
    opex = length(producers) * mdp.opex

    return revenue - capex - opex
end


function POMDPs.initialstate(mdp::GeothermalMDP)
    return Deterministic(State(t = mdp.t0, i = [], p = []))
end

function POMDPs.discount(mdp::GeothermalMDP)
    return mdp.γ
end