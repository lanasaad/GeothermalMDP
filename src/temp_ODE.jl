function temperature_model(df::DataFrame, injectors::Vector{Coord}, producers::Vector{Coord}, simulation_time::Float64, flow_rate::Float64=8000.0, t_in::Float64=30.0, ambient_radius::Float64=30.0, ambient_temp::Float64=60.0)
    function temp_ode!(du, u, p, t)
        df_temp = u
        for i in 1:length(df_temp)

            #TODO: calculate ambient temperature dynamically
            # ambient_temps = []
            # for j in 1:length(df_temp)
            #     distance = euclidean([df.x[i], df.y[i]], [df.x[j], df.y[j]])
            #     if distance <= ambient_radius
            #         push!(ambient_temps, df_temp[j])
            #     end
            # end
            
            # # If no nearby temperatures, use the original 20°C
            # ambient_temp = length(ambient_temps) > 0 ? mean(ambient_temps) : ambient_temp

            # Natural dissipation
            du[i] = -0.1 * (df_temp[i] - ambient_temp)  
            
            # Effect of injectors
            for inj in injectors
                distance = euclidean([df.x[i], df.y[i]], [inj.x, inj.y])
                du[i] += -0.01 * flow_rate * (df_temp[i] - t_in) / (distance^2 + 1)
                # println("Inj Distance: ", distance)
            end
            
            # Effect of producers
            for prod in producers
                distance = euclidean([df.x[i], df.y[i]], [prod.x, prod.y])
                du[i] += -0.01 * flow_rate * (df_temp[i] - df_temp[i]) / (distance^2 + 1)
                # println("Inj Distance: ", distance)
            end
        end
    end

    tspan = (0.0, simulation_time)
    prob = ODEProblem(temp_ode!, df.temperature, tspan)
    sol = DifferentialEquations.solve(prob, Tsit5())
    
    return sol[end]
end


s = State(t = df.temperature, i=[mdp.well_pairs[1][1]], p=[mdp.well_pairs[1][2]])
flow_rate = 8000.0
t_in = 30.0
simulation_time = 10.0 
df_next = temperature_model(df, s.i, s.p, simulation_time)

plot_temp(df, df.temperature)
plot_temp(df, df_next)

