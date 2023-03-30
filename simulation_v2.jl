using Dates
using DifferentialEquations
using Sundials
# using ODEInterfaceDiffEq
using Plots
using BenchmarkTools
using Statistics
using SparseArrays
using LinearAlgebra

# using ModelingToolkit


using MasterChemicalMechanism

fpath = "./src/data/extracted/alkanes/methane.fac"
# fpath = "./src/data/extracted/alkanes/meth_eth_prop_but.fac"
# fpath = "./src/data/extracted/monoterpines/alpha_pinene.fac"
# fpath = "./src/data/extracted/monoterpines/limonene.fac"
# fpath = "./src/data/extracted/monoterpines/monoterpines.fac"
# fpath = "./src/data/extracted/alkanes/all_alkanes.fac"
# fpath = "./src/data/extracted/full/mcm_subset.fac"
# fpath = "./src/data/extracted/no_terpenes.fac"
# fpath = "./src/data/extracted/through_but.fac"

model_name = "methane"
@assert ispath(fpath)  == true

fac_dict = read_fac_file(fpath)


keys(fac_dict)
fac_dict["reaction_definitions"][1]


generate_config(fac_dict; model_name=model_name)
include("./model/$(model_name)/config.jl")  # <-- we need to update to automatically generate this


rate_list = generate_rrates(fac_dict; model_name=model_name, params=starting_params)
include("./model/$(model_name)/rrates.jl")


generate_species(fac_dict; model_name=model_name)
include("./model/$(model_name)/species.jl")


generate_ro2(fac_dict; model_name=model_name)
include("./model/$(model_name)/ro2.jl")


generate_photolysis_mcm(fac_dict; model_name)
include("./model/$(model_name)/photolysis.jl")


generate_stoich_mat(fac_dict; model_name=model_name)
# overall, reactants, reactant_indicator

const N = get_sparse_mat(mat_name="N", model_name=model_name)
const Nrows = rowvals(N)
const Nvals = nonzeros(N)

const R = get_sparse_mat(mat_name="R", model_name=model_name)
const Rrows = rowvals(R)
const Rvals = nonzeros(R)


generate_rrates_mechanism(fac_dict, rate_list; model_name=model_name, params=starting_params)
include("./model/$(model_name)/rrates_mechanism.jl")
const k_rates = zeros(size(R,2))
const v = zeros(size(R,2))
const Jtemp = 0.0

size(R)
idxs_reactants = findall(>(0), R[:,1])

# compute conversion to ppb from molecules/cc
m_init = M(starting_params.T, starting_params.P)
nₘ = 10^9 # for ppb

# generate sane initial conditions
init_path = "./src/data/initial_concentrations/full.txt"
isfile(init_path)
init_dict = generate_init_dict(init_path, m_init)

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
const RO2 = sum(u₀[idx_ro2])

n_days = 1
tspan = (0.0, n_days*24.0*60.0)
#tspan = (0.0, 60.0)
#tspan = (0.0, 15.0)
tol = 1e-6


# -----------------------------------------
# generate derivative terms
# -----------------------------------------
const derivative_terms = get_derivative_terms(N)
const reaction_terms = get_reaction_terms(R)

# -----------------------------------------
# copute jacobian terms
# -----------------------------------------

const Jprototype = get_jac_prototype(N,R)
const jac_terms = get_jac_terms(Jprototype, N, R)

size(derivative_terms)
size(reaction_terms)
size(jac_terms)

jac_terms[1]

reaction_terms[1]

const du_temp = 0.0
const Jtemp = 0.0

# ps = (starting_params.T, starting_params.P)

const T = starting_params.T
const P = starting_params.P
# ps = (starting_params.T, starting_params.P, N, R, derivative_terms, reaction_terms, idx_ro2, RO2, k_rates, du_temp, jac_terms, Jtemp)
# ps = (starting_params.T, starting_params.P, N, R, Rrows, Rvals, idx_ro2, RO2, k_rates, v)


# take a sample step to make sure everything is pre-compiled all nice
u₀
du = copy(u₀)
# f!(du, u₀, ps, 0.0)
f!(du, u₀, nothing, 0.0)
du
@benchmark f!(du, u₀, nothing, 0.0)



# Jtest = copy(Jprototype)
# Jac!(Jtest, u₀, ps, 0.0)
# Jtest

# @benchmark Jac!(Jtest, u₀, ps, 0.0)


# set up second odeproblem for comparison
#fun = ODEFunction(f!; jac=Jac!, jac_prototype=JP)
ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(f!, u₀, tspan)
#ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(f!, u₀, tspan, ps)
# fun = ODEFunction(f!; jac=Jac!, jac_prototype=Jprototype)
# fun = ODEFunction(f!;  jac_prototype=Jprototype)
# ode_prob_2 = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun, u₀, tspan, ps)


using Zygote, SciMLSensitivity


# we should test forward sensitivity and reverse via QuadratureAdjoint
# ForwardSenitivity
# QuadratureAdjoint
# BacksolveAdjoint
# InterpolatingAdjoint

function sum_of_solution(u₀)
    _prob = remake(ode_prob; u0=u₀, p=[0.0])
    sum(solve(_prob, CVODE_BDF(); saveat=15.0, reltol=tol, abstol=tol, sensealg=QuadratureAdjoint()))
end

function sum_of_solution2(u₀)
    _prob = remake(ode_prob; u0=u₀, p=[0.0])
    sum(solve(_prob, CVODE_BDF(); saveat=15.0, reltol=tol, abstol=tol, sensealg=InterpolatingAdjoint()))
end


function sum_of_solution3(u₀)
    _prob = remake(ode_prob; u0=u₀, p=[0.0])
    sum(solve(_prob, CVODE_BDF(); saveat=15.0, reltol=tol, abstol=tol, sensealg=BacksolveAdjoint()))
end

@benchmark Zygote.gradient(sum_of_solution, u₀)
@benchmark Zygote.gradient(sum_of_solution2, u₀)
# @benchmark Zygote.gradient(sum_of_solution3, u₀) unstable




@benchmark solve(
    ode_prob,
    CVODE_BDF();
    saveat=15.0,
    reltol=tol,
    abstol=tol,
)

@benchmark solve(
    ode_prob_2,
    CVODE_BDF();
    saveat=15.0,
    reltol=tol,
    abstol=tol,
)

# @benchmark solve(
#     ode_prob_2,
#     QNDF();
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

