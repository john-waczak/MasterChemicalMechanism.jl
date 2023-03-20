using ModelingToolkit
using Catalyst
using Dates
using Latexify
using DifferentialEquations

using MasterChemicalMechanism
import MasterChemicalMechanism: read_fac_file, parse_rxns
import MasterChemicalMechanism: generate_config
import MasterChemicalMechanism: generate_rrates
import MasterChemicalMechanism: generate_species
import MasterChemicalMechanism: generate_ro2
import MasterChemicalMechanism: generate_photolysis_mcm
import MasterChemicalMechanism: generate_mechanism

using Plots
using BenchmarkTools


# fpath = "./src/data/extracted/alkanes/methane.fac"
# fpath = "./src/data/extracted/alkanes/meth_eth_prop_but.fac"
# fpath = "./src/data/extracted/monoterpines/alpha_pinene.fac"
#fpath = "./src/data/extracted/monoterpines/limonene.fac"
# fpath = "./src/data/extracted/monoterpines/monoterpines.fac"
# fpath = "./src/data/extracted/alkanes/all_alkanes.fac"
fpath = "./src/data/extracted/full/mcm_subset.fac"
# fpath = "./src/data/extracted/no_terpenes.fac"
#fpath = "./src/data/extracted/through_but.fac"
model_name = "mcm_full"
@assert ispath(fpath)  == true


# 1. constants
# 2. reaction rates
fac_dict = read_fac_file(fpath)
#rxns = fac_dict["reaction_definitions"]
#species, reactions = parse_rxns(rxns)

generate_config(fac_dict, model_name)
generate_rrates(fac_dict, model_name)
generate_species(fac_dict, model_name)
generate_ro2(fac_dict, model_name)
generate_photolysis_mcm(fac_dict, model_name)
generate_mechanism(fac_dict, model_name)

include("./model/$(model_name)/config.jl")  # <-- we need to update to automatically generate this
include("./model/$(model_name)/rrates.jl")
include("./model/$(model_name)/species.jl")
include("./model/$(model_name)/ro2.jl")

include("./model/$(model_name)/photolysis.jl")
# test that the photolysis works
times = 0.0:1.0:(60*60*24)
j1 = J_1.(times)

include("./model/$(model_name)/mechanism.jl")


# ----------------------------------------------------------------
# Running the simulation
# ----------------------------------------------------------------

#ndays = 2
tspan = (0.0, 60.0)
# use these initial concentrations and set everything
# else to 0.0: https://github.com/AtChem/AtChem2/blob/master/model/configuration/initialConcentrations.config

# init_dict = Dict(
#     "CH4" => 4.9e+13,
#     "CO" => 3.6e+12,
#     "O3" => 5.2e+11,
#     "NO2" => 2.4e+11,
# )

m_init = substitute(M, Dict([
    T => 291.483, # 65 °F
    P => 1013.2,  # milibars standard pressure
    kb => 1.380649e−23 # J/K
]))

m_init = m_init.val
nₘ = 10^9 # for ppb

# convert ppb to molecules/cc
init_dict = Dict(
    "CH4" => 1909.6*m_init/nₘ,
    "CO" => 100.0*m_init/nₘ,
    "O3" => 18.0*m_init/nₘ,
    "NO2" => 100*m_init/nₘ,
    "APINENE" => 30.0*m_init/nₘ,
    "H2" => 500*m_init/nₘ,
)


# see this link: https://docs.sciml.ai/Catalyst/stable/example_networks/smoluchowski_coagulation_equation/#smoluchowski_coagulation_equation
u₀    = zeros(Float64, size(X,1))
for (key, val) ∈ init_dict
    try
        println("$(key): $(val)")
        idx = findfirst(x -> x == key, species)
        println("idx: ", idx)
        u₀[idx] = val
    catch e
        println(e)
    end
end
u₀map = Pair.(collect(X), u₀)   # map variable to its initial value

params = (
    T => 291.483, # 65 °F
    P => 1013.2  # milibars standard pressure
)


# let's convert the model into an ODESystem so we can simplify

odesys = @time convert(ODESystem, mechanism; combinatoric_ratelaws=false)
ode_prob = @time ODEProblem(odesys, u₀, tspan, params)  # 600.706338 seconds


fieldnames(typeof(ode_prob))

du = copy(u₀)
ode_prob.f(du, u₀, (291.483, 1013.2), 0.0)


# sol = solve(ode_prob, saveat=15);
println(fpath)
println("num species: ", numspecies(mechanism))  # 3490
println("num reactions: ", numreactions(mechanism)) # 1163

# @btime solve(ode_prob;
#             alg_hints=[:stiff],
#             dense=false,
#             saveat=15,
#             reltol=1e-6,
#             )

# solve(ode_prob;
#       alg_hints=[:stiff],
#       dense=false,
#       saveat=15.0,
#       save_on=false,  # nuclear option
#       reltol=1e-6,
#       )

integrator = init(ode_prob,
                  TRBDF2();
#                  alg_hints=[:stiff],
                  dense=false,
                  save_everystep=false,
                  saveat=15.0,
#                  save_on=false,  # nuclear option
                  reltol=1e-6,

                  )

step!(integrator)

varinfo()


#using LSODA
#sol = @time solve(ode_prob, lsoda(), saveat=15)
#sol = solve(ode_prob, AutoTsit5(Rosenbrock23()), saveat=15)  # <-- didn't work
#sol = solve(ode_prob, TRBDF2(), saveat=15)
#sol = solve(ode_prob, QNDF(), saveat=15)
#sol = solve(ode_prob, FBDF(), saveat=15)



# unique_species

# speciesiwant = ["CO", "O3", "NO2", "APINENE", "H2"]
# idxs= [findfirst(x->x==spec, unique_species) for spec ∈ speciesiwant]

# p = plot()
# i = 1
# for idx ∈ idxs
#     try
#         plot!(p, sol.t ./ (60*24), sol[idx, :] .* (nₘ/m_init), label=speciesiwant[i], xlabel="t [days]", ylabel="concentration [ppb]", legend=:outertopright)
#     catch e
#         println(i, idxs[i])
#     end

#     i += 1
# end


# display(p)
