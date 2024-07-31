using POMDPs  
using POMDPModelTools
using POMDPPolicies
using Parameters
using Random
using CSV
using DataFrames
using Revise
using IterTools


#won't use now as we're looking at it at small scale
#file_path = "/Users/lanas/Desktop/GeothermalPOMDP/res_temppress_t5_k184.csv"
#data = CSV.read(file_path, DataFrame)

@with_kw mutable struct State
    temperatures::Array{Float64,2} 
    well_states::Array{Int,2}
end

 
@with_kw mutable struct GeoCell
    i::Int
    j::Int
end


@with_kw mutable struct DrillAction
    injector::GeoCell
    producer::GeoCell
end


@with_kw mutable struct ShutAction
    injector::GeoCell
    producer::GeoCell
end


@with_kw mutable struct NothingAction
    injector::GeoCell=GeoCell(-1, -1)
    producer::GeoCell=GeoCell(-1, -1)
end


const Action = Union{DrillAction,ShutAction, NothingAction}

const NO_WELL = 0
const INJECTION_WELL = 1
const PRODUCTION_WELL = 2



#initial States
well_states_grid = fill(NO_WELL, (5, 5))
temperatures_grid::Array{Float64, 2} = fill(100.0, 5, 5)



temperatures_list=[]
is_well_list=[]

for _ in 1:100
    push!(temperatures_list, rand([40,60.0,80.0,100.0,120.0],5, 5))
    push!(is_well_list,rand([NO_WELL,INJECTION_WELL,PRODUCTION_WELL],5,5))
end

print(temperatures_list)
temperatures_list[1]
is_well_list[1]
#typeof(temperature_list)



#specify injector and producer locations here to populate state space
@with_kw mutable struct GeothermalMDP <: MDP{State , Action}
	size::Tuple{Int,Int} = (5,5)   # size of the grid
    initial_temperatures::Array{Float64, 2} = temperatures_grid  
    initial_well_states::Array{Int, 2} = well_states_grid
    temp_states::Vector{Any} = temperatures_list
    is_well_states::Vector{Any} = is_well_list

    injector_positions::Vector{Tuple{Int, Int}} = [(1, 1), (2, 2),(3,1),(4,2),(5,1)]  
    producer_positions::Vector{Tuple{Int, Int}} = [(1, 3), (2, 4),(3,3),(4,4),(5,5)] 

    Δtemp::Float64= 20.0
    
    #defining values for capex and opex?
    capex::Float64 = 10000 #for each well pair
    opex::Float64 = 1000 #for each well pair
    revenue::Float64 = 5000

    γ::Real = 0.95
end	



#states includes the coordinates in the form of x,y and temperatures across the field 
function POMDPs.states(mdp::GeothermalMDP)
    
    #[(i, j, mdp.temperatures[i, j], mdp.well_states[i, j]) for i in 1:mdp.size[1], j in 1:mdp.size[2]]
    states=[]

    for i in 1:length(mdp.temp_states)
        s = State(temperatures=mdp.temp_states[i],well_states_grid=mdp.is_well_states[i])
        push!(states,s)
    end 
    return states


    #TODO fiz states
    range_temp = [40,60.0,80.0,100.0,120.0]
    range_well_state = [0,1,2]
    𝒮 = [State(s[1], s[2]) for s in Iterators.product(range_temp, range_well_state)]
    return 𝒮

end


#states(g)



#Actions for mdp
function POMDPs.actions(mdp::GeothermalMDP)

   injector_locations=[[1,2],[3,4],[2,5]]
   producer_locations=[[1,1],[2,3],[4,4]]

   actions=[]
   

   for k in 1:3
    d=DrillAction(
        injector=GeoCell(i=injector_locations[k][1],j=injector_locations[k][2]),
        producer=GeoCell(i=producer_locations[k][1],j=producer_locations[k][2])
        
        )
    s=ShutAction(
            injector=GeoCell(i=injector_locations[k][1],j=injector_locations[k][2]),
            producer=GeoCell(i=producer_locations[k][1],j=producer_locations[k][2])
            
            )

    push!(actions,d)
    push!(actions,s)
   end

   nothing_action = NothingAction()

   push!(actions, nothing_action)
   return actions
end




#TODO assumption Drill->Operate

function POMDPs.transition(mdp::GeothermalMDP,s::State, a::Action)
    sp=deepcopy(s) 
    #PLANNING step
    println("Initial copied state: ", sp)  # Print the initial state copy

    if typeof(a)==DrillAction
        inj_i=a.injector.i
        inj_j=a.injector.j
        prod_i=a.producer.i
        prod_j=a.producer.j

        sp.well_states[inj_i,inj_j]=INJECTION_WELL
        sp.well_states[prod_i,prod_j]=PRODUCTION_WELL


    elseif typeof(a)==ShutAction
            inj_i=a.injector.i
            inj_j=a.injector.j
            prod_i=a.producer.i
            prod_j=a.producer.j
    
    
            sp.well_states[inj_i,inj_j]=NO_WELL
            sp.well_states[prod_i,prod_j]=NO_WELL


#operation keeps going for existing wells 
    end 

    #operations
    injectors = findall(x -> x==INJECTION_WELL, sp.well_states)
    producers = findall(x -> x==PRODUCTION_WELL,sp.well_states)

    all_wells=vcat(injectors,producers)


    wtemp_list=[]
    wtemp_probs=[]
    wttemp_template = 
    Δtemps = [0, -mdp.Δtemp, -2 * mdp.Δtemp]
    Δtemp_probs = [0.7, 0.25, 0.05] #TODO: move to the problem struct
    

    for well_inj in injectors
        for well_prod in producers
            for (i, Δtemp) in enumerate(Δtemps)
                new_temps = deepcopy(sp.temperatures)
                new_temps[well_inj.I[1], well_inj.I[2]] += Δtemp
                new_temps[well_prod.I[1], well_prod.I[2]] += Δtemp
                push!(wtemp_list, new_temps)
            end
        end
    end


        
        
        new_temperatures= [current_temperatures, current_temperatures - mdp.Δtemp, current_temperatures - 2 * mdp.Δtemp]
        

        #TODO create a SparseCat for all combinations of new temperature values (prod)
        temperature_distribution = SparseCat(new_temperatures, probabilities)

        new_temperatures=Dict("i"=>well.I[1],
        "j"=>well.I[2],
        "values"=>new_temperatures)
        push!(well_temperatures,new_temperatures)


        #TODO returning rand temp out of this distribution


        #sp.temperatures[well.I[1],well.I[2]]= mdp.Δtemp

    
    end 



    return Deterministic(sp)
end




#Operation step- for all existing wells, update temperature
#i, j = injectors[k].I


function POMDPs.reward(mdp::GeothermalMDP,s::State,a::Action,sp::State)
    #rev * no.of prod

    #operations
    injectors = findall(x -> x==INJECTION_WELL, sp.well_state)
    producers = findall(x -> x==PRODUCTION_WELL,sp.well_state)

    num_injectors=length(injectors)
    num_producers=length(producers)

    revenue=num_producers*mdp.revenue
    capex=typeof(a)==DrillAction ? mdp.capex : 0
    opex=num_injectors*mdp.opex

    r = revenue-capex-opex

    return r
end 


function POMDPs.initialstate(mdp::GeothermalMDP)
    return State(temperatures=mdp.initial_temperatures,well_states=mdp.initial_well_states)
end












































#=
# Assuming you have defined your GeothermalMDP and other necessary components
b = GeothermalMDP()

# Get the initial state
initial_state = POMDPs.initialstate(b)
# To display or use the initial state, you can print it out or inspect its properties


=#

