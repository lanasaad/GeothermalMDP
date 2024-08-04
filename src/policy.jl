struct BuildRandom <: Policy
	g::GeothermalMDP
end

function POMDPs.action(p::BuildRandom, s::State)
	if length(s.i) < 1
		return rand(rng, actions(p.g)[1:1000])
	else
		return rand(rng, actions(p.g, s))
	end
end
