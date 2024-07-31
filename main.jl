using POMDPSimulators
using POMDPTools
using POMDPPolicies

include("pomdp.jl")

g=GeothermalMDP()

# rand_policy(g::GeothermalMDP, s::State) = rand(actions(g))

#first iteration
s = initialstate(g)
#a = rand_policy(g, s)
policy = RandomPolicy(g)

actions(g)
a = action(policy, s)
a.injector, a.producer

sp = rand(transition(g,s,a))
r = reward(g,s,a,sp)
sum(sp.well_states)

injectors = findall(x -> x==INJECTION_WELL, sp.well_states)
producers = findall(x -> x==PRODUCTION_WELL,sp.well_states)

w1 = injectors[1]
typeof(w1)
#2nd iteration
s = deepcopy(sp)
a = rand_policy(g, s)
sp = rand(transition(g,s,a))
r = reward(g,s,a,sp)


s = deepcopy(sp)
a = rand_policy(g, s)
sp = rand(transition(g,s,a))
r = reward(g,s,a,sp)

c = states(g)