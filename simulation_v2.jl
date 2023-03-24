using Dates
using DifferentialEquations
using Sundials
# using ODEInterfaceDiffEq
using Plots
using BenchmarkTools
using Statistics
using SparseArrays
using LinearAlgebra

using MasterChemicalMechanism

# fpath = "./src/data/extracted/alkanes/methane.fac"
# fpath = "./src/data/extracted/alkanes/meth_eth_prop_but.fac"
# fpath = "./src/data/extracted/monoterpines/alpha_pinene.fac"
fpath = "./src/data/extracted/monoterpines/limonene.fac"
# fpath = "./src/data/extracted/monoterpines/monoterpines.fac"
# fpath = "./src/data/extracted/alkanes/all_alkanes.fac"
# fpath = "./src/data/extracted/full/mcm_subset.fac"
# fpath = "./src/data/extracted/no_terpenes.fac"
# fpath = "./src/data/extracted/through_but.fac"

model_name = "limonene"
@assert ispath(fpath)  == true

fac_dict = read_fac_file(fpath)


keys(fac_dict)
fac_dict["reaction_definitions"][1]


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


generate_stoich_mat(fac_dict; model_name=model_name)
# overall, reactants, reactant_indicator
N = get_sparse_mat(mat_name="N", model_name=model_name)
R = get_sparse_mat(mat_name="R", model_name=model_name)
R_idxs = [findall(>(0), R[:,j]) for j ∈ axes(R,2)]


generate_rrates_mechanism(fac_dict, rate_list; model_name=model_name, params=params)
include("./model/$(model_name)/rrates_mechanism.jl")
k_rates = zeros(size(R,2))
v = zeros(size(R,2))
Jtemp = 0.0

size(R)
idxs_reactants = findall(>(0), R[:,1])

# generate_ode_f(fac_dict, idx_ro2; model_name=model_name, params=params)
# include("./model/$(model_name)/ode_f.jl")
# # include("./model/$(model_name)/ode_jac.jl")


# compute conversion to ppb from molecules/cc
m_init = M(params.T, params.P)
nₘ = 10^9 # for ppb

# generate sane initial conditions
init_path = "./src/data/initial_concentrations/full.txt"
isfile(init_path)
init_dict = generate_init_dict(init_path, m_init)


# # convert ppb to molecules/cc
# init_dict = Dict(
#     "CH4" => 1909.6*m_init/nₘ,
#     "CO" => 100.0*m_init/nₘ,
#     "O3" => 18.0*m_init/nₘ,
#     "NO2" => 100*m_init/nₘ,
#     "APINENE" => 30.0*m_init/nₘ,
#     "H2" => 500*m_init/nₘ,
# )


# see this link: https://docs.sciml.ai/Catalyst/stable/example_networks/smoluchowski_coagulation_equation/#smoluchowski_coagulation_equation
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

# combine parameters into one long tuple
RO2 = sum(u₀[idx_ro2])
ps = (params.T, params.P, N, R, R_idxs, idx_ro2, RO2, k_rates, v, Jtemp)


n_days = 1
tspan = (0.0, n_days*24.0*60.0)
#tspan = (0.0, 15.0)
tol = 1e-6
#ode_prob = @time ODEProblem(mcm!, u₀, tspan, ps)  # 600.706338 seconds
ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(f!, u₀, tspan, ps)  # 600.706338 seconds

# take a sample step to make sure everything is pre-compiled all nice
u₀
du = copy(u₀)
@btime f!(du, u₀, ps, 0.0)
f!(du, u₀, ps, 0.0)

Jtest = zeros(size(R,1), size(R,1))
@btime J!(Jtest, u₀, ps, 0.0)

# now that we've evaluated the jacobian, generate a sparse prototype to use:
J!(Jtest, u₀, ps, 0.0)
JP = sparse(Jtest)
JP .= 0.0

# set up second odeproblem for comparison
fun = ODEFunction(f!; jac=J!, jac_prototype=JP)
ode_prob_2 = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun, u₀, tspan, ps)



# @benchmark solve(
#     ode_prob,
#     CVODE_BDF();
#     saveat=15.0,
#     reltol=tol,
#     abstol=tol,
# )

# @benchmark solve(
#     ode_prob_2,
#     CVODE_BDF();
#     saveat=15.0,
#     reltol=tol,
#     abstol=tol,
# )




sol = @time solve(
    ode_prob,
    CVODE_BDF();
    saveat=15.0,
    reltol=tol,
    abstol=tol,
);

sol2 = @time solve(
    ode_prob_2,
    CVODE_BDF();
    saveat=15.0,
    reltol=tol,
    abstol=tol,
);




# # calculate mean dt
# mean_dt = mean(sol.t[2:end] .- sol.t[1:end-1])
# std_dt = std(sol.t[2:end] .- sol.t[1:end-1])


# solver_name = "CVODE_BDF"
# #solver_name = "ODEX"
# p = plot(sol.t[2:end] .- sol.t[1:end-1], lw=2, xlabel="step", ylabel="dt [minutes]", title=solver_name, label="")
# savefig("dt_$(solver_name)_$(n_days)-days_nojac.pdf")
# savefig("dt_$(solver_name)_$(n_days)-days_nojac.png")
# display(p)




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


p2 = plot()
i = 1
for idx ∈ idxs
    try
        plot!(p2, sol2.t ./ (60*24), sol2[idx, :] .* (nₘ/m_init), label=speciesiwant[i], xlabel="t [days]", ylabel="concentration [ppb]", legend=:outertopright)
    catch e
        println(i, idxs[i])
    end

    i += 1
end

#dsipaly(p)

plot(p,p2)

ylims!(0.0, 100.0)

