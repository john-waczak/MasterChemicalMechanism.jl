module MasterChemicalMechanism

using ModelingToolkit
using Catalyst
using SolarGeometry

include("parse_fac.jl")
include("config.jl")
include("reaction_rates.jl")
include("species.jl")
include("ro2.jl")
include("photolysis.jl")
include("mechanism.jl")

# 6. define photolysis rates
# include("photolysis.jl")

# 7. generate reaction_network via catalyst
# include("mechanism.jl")


# include("catalyst_model.jl")
# include("simulate.jl")


#we want to create a function that takes the .fac file and spits out a julia file defining the entire mechanism via catalyst, modelingtoolkit. If needed, we can copy the files


end



# use Event Handling and Callback Functions to update parameters (e.g. T, P, RH, etc...) with measurement values
# also use it to update the state vector (i.e. concentrations, reaction rate functions, photolysis rates, ROx)
# via the data assimilation step for Extended Kalman Filter
# https://docs.sciml.ai/ModelingToolkit/stable/basics/Events/




# m=2.55e+19
# o2=0.2095*m


# REAL(dp)::M, N2, O2, RO2, H2O, TEMP,

# We should specify
# temp, pressure, H2O as variables which we can set initially or read from datafile


# RO2 list: https://github.com/AtChem/AtChem2/blob/master/mcm/peroxy-radicals_v3.3.1
# photolysis rates: https://github.com/AtChem/AtChem2/blob/master/mcm/photolysis-rates_v3.3.1

# atmospheric functions (i.e. M calculation) https://github.com/AtChem/AtChem2/blob/37503aefb02bb54fa3b04899c68d15afa8b1c8c0/src/atmosphereFunctions.f90
