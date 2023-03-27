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
N = get_sparse_mat(mat_name="N", model_name=model_name)
Nrows = rowvals(N)
Nvals = nonzeros(N)

R = get_sparse_mat(mat_name="R", model_name=model_name)
Rrows = rowvals(R)
Rvals = nonzeros(R)
Vmat = copy(float(R)') # we should make sure this isn't a union type


generate_rrates_mechanism(fac_dict, rate_list; model_name=model_name, params=starting_params)
include("./model/$(model_name)/rrates_mechanism.jl")
k_rates = zeros(size(R,2))
v = zeros(size(R,2))
Jtemp = 0.0

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
RO2 = sum(u₀[idx_ro2])

n_days = 1
tspan = (0.0, n_days*24.0*60.0)
#tspan = (0.0, 15.0)
tol = 1e-6



# -----------------------------------------
# generate derivative terms
# -----------------------------------------
derivative_terms = get_derivative_terms(N)
reaction_terms = get_reaction_terms(R)

derivative_terms
reaction_terms
# -----------------------------------------
# copute jacobian terms
# -----------------------------------------

jac_terms = JacobianTerms[]
Jprototype = zeros(Float64, size(R,1), size(R,1))

In,Kn,Vn = findnz(N)

for (i,k) ∈ zip(In, Kn)
    for j ∈ axes(R,1)
        if N[i,k] != 0 && R[j,k] > 0
            # isinJ[i,j,k] = 1
            Jprototype[i,j] = 1.0
        end
    end
end


# Jprototype2 = zeros(Float64, size(R,1), size(R,1))
# for k ∈ axes(R,2)
#     for j ∈ axes(R,1)
#         for i ∈ axes(N,1)
#             if N[i,k] != 0 && R[j,k] > 0
#                 # isinJ[i,j,k] = 1
#                 Jprototype2[i,j] = 1.0
#             end
#         end
#     end
# end

# all(Jprototype .== Jprototype2)


Jprototype = sparse(Jprototype)
println("nonzero %: ", 100*length(nonzeros(Jprototype))/(size(Jprototype,1)^2))

# loop over all (i,j) pairs in Jprototype

Jprot_rows = rowvals(Jprototype)
Jprot_vals = nonzeros(Jprototype)

for j ∈ axes(Jprototype, 1)
    for i_row ∈ nzrange(Jprototype, j)
        i = Jprot_rows[i_row]
        rxn_terms = ReactionTerms[]
        for k ∈ axes(R,2)
            if N[i,k] != 0 && R[j,k] > 0
                reaction_indices = findall(x -> x > 0, R[:,k])
                @assert j ∈ reaction_indices
                indices_out = [j]
                # add the rest of the indices so it's sorted with j first
                for idx_rxn ∈ reaction_indices
                    if idx_rxn != j
                        push!(indices_out, idx_rxn)
                    end
                end

                push!(rxn_terms, ReactionTerms(k, indices_out))
            end
        end
        push!(jac_terms, JacobianTerms(i,j,rxn_terms))
    end
end

size(derivative_terms)
size(reaction_terms)
size(jac_terms)

jac_terms[1]

reaction_terms[1]

du_temp = 0.0
Jtemp = 0.0
ps = (starting_params.T, starting_params.P, N, R, derivative_terms, reaction_terms, idx_ro2, RO2, k_rates, du_temp, jac_terms, Jtemp)

# take a sample step to make sure everything is pre-compiled all nice
u₀
du = copy(u₀)
f!(du, u₀, ps, 0.0)
du
@benchmark f!(du, u₀, ps, 0.0)

Jtest = copy(Jprototype)
Jac!(Jtest, u₀, ps, 0.0)
Jtest

@benchmark Jac!(Jtest, u₀, ps, 0.0)



# set up second odeproblem for comparison
#fun = ODEFunction(f!; jac=Jac!, jac_prototype=JP)
ode_prob = @time ODEProblem{true, SciMLBase.FullSpecialize}(f!, u₀, tspan, ps)
fun = ODEFunction(f!; jac=Jac!, jac_prototype=Jprototype)
ode_prob_2 = @time ODEProblem{true, SciMLBase.FullSpecialize}(fun, u₀, tspan, ps)

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

