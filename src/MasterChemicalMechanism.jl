module MasterChemicalMechanism

include("atmospheric_functions.jl")
# include("config.jl")
# include("mechanism.jl")
# include("reaction_rates.jl")
# include("catalyst_model.jl")
# include("simulate.jl")
# include("photolysis.jl")

end



# use Event Handling and Callback Functions to update parameters (e.g. T, P, RH, etc...) with measurement values
# also use it to update the state vector (i.e. concentrations, reaction rate functions, photolysis rates, ROx)
# via the data assimilation step for Extended Kalman Filter
# https://docs.sciml.ai/ModelingToolkit/stable/basics/Events/
