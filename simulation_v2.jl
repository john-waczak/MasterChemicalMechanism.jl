using Dates
using DifferentialEquations
using Sundials
using Plots
# using BenchmarkTools
#using MasterChemicalMechanism

import MasterChemicalMechanism: read_fac_file, parse_rxns
import MasterChemicalMechanism: generate_config
import MasterChemicalMechanism: generate_rrates
import MasterChemicalMechanism: generate_species
import MasterChemicalMechanism: generate_ro2
import MasterChemicalMechanism: generate_photolysis_mcm
import MasterChemicalMechanism: generate_rrates_mechanism
import MasterChemicalMechanism: generate_ode_f


# fpath = "./src/data/extracted/alkanes/methane.fac"
# fpath = "./src/data/extracted/alkanes/meth_eth_prop_but.fac"
# fpath = "./src/data/extracted/monoterpines/alpha_pinene.fac"
# fpath = "./src/data/extracted/monoterpines/limonene.fac"
# fpath = "./src/data/extracted/monoterpines/monoterpines.fac"
# fpath = "./src/data/extracted/alkanes/all_alkanes.fac"
fpath = "./src/data/extracted/full/mcm_subset.fac"
# fpath = "./src/data/extracted/no_terpenes.fac"
# fpath = "./src/data/extracted/through_but.fac"
model_name = "mcm_full"
@assert ispath(fpath)  == true

fac_dict = read_fac_file(fpath)



generate_config(fac_dict; model_name=model_name)
include("./model/$(model_name)/config.jl")  # <-- we need to update to automatically generate this


rate_list = generate_rrates(fac_dict; model_name=model_name, params=params)
include("./model/$(model_name)/rrates.jl")


generate_species(fac_dict; model_name=model_name)
include("./model/$(model_name)/species.jl")


generate_ro2(fac_dict; model_name=model_name)
include("./model/$(model_name)/ro2.jl")


generate_photolysis_mcm(fac_dict; model_name)
include("./model/$(model_name)/photolysis.jl")


generate_rrates_mechanism(fac_dict, rate_list; model_name=model_name, params=params)
include("./model/$(model_name)/rrates_mechanism.jl")


generate_ode_f(fac_dict, idx_ro2; model_name=model_name, params=params)
include("./model/$(model_name)/ode_f.jl")



# combine parameters into one long tuple
ps = (params.T, params.P, idx_ro2)


# compute conversion to ppb from molecules/cc
m_init = M(params.T, params.P)
nₘ = 10^9 # for ppb


# convert ppb to molecules/cc
init_dict = Dict(
    "CH4" => 1909.6*m_init/nₘ,
    "CO" => 100.0*m_init/nₘ,
    "O3" => 18.0*m_init/nₘ,
    "NO2" => 100*m_init/nₘ,
#    "APINENE" => 30.0*m_init/nₘ,
    "H2" => 500*m_init/nₘ,
)


# see this link: https://docs.sciml.ai/Catalyst/stable/example_networks/smoluchowski_coagulation_equation/#smoluchowski_coagulation_equation
#u₀    = (1.0*m_init/nₘ) .* ones(Float64, size(species))
u₀    = zeros(Float64, size(species))
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

n_days = 5
tspan = (0.0, n_days*24.0*60.0)

ode_prob = @time ODEProblem(mcm!, u₀, tspan, ps)  # 600.706338 seconds

# u₀
# du = copy(u₀)
# mcm!(du, u₀, ps, 0.0)



# sol = solve(ode_prob, saveat=15);
# @btime solve(ode_prob;
#             alg_hints=[:stiff],
#             dense=false,
#             saveat=15,
#             reltol=1e-6,
#             )

sol = @time solve(ode_prob,
                  CVODE_BDF();
                  dense=false,
                  saveat=15.0,
                  reltol=1e-5,
                  abstol=1e-5,
#                  reltol=1e-6,
#                  abstol=1e-6,
                  )


# integrator = init(ode_prob;
#                   alg_hints=[:stiff],
#                   dt=dt,
#                   saveat=15.0,
#                   dense=false,
#                   reltol=1e-6
#                   )

# @time step!(integrator)



println(size(species))

#speciesiwant = species[1:10]

speciesiwant = species[1:20]

#speciesiwant = ["CO", "O3", "NO2", "APINENE", "H2", "CH4"]
#speciesiwant = ["H2O2"]

idxs= [findfirst(x->x==spec, species) for spec ∈ speciesiwant]


p = plot()
i = 1
for idx ∈ idxs
    try
        plot!(p, sol.t ./ (60*24), sol[idx, :] .* (nₘ/m_init), label=speciesiwant[i], xlabel="t [days]", ylabel="concentration [ppb]", legend=:outertopright)
    catch e
        println(i, idxs[i])
    end

    i += 1
end

# ylims!(0.0, 100.0)
# ylims!(0.0, 100.0)

display(p)
