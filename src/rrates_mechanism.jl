function rxnrate_to_expression(rxnrate_string, rate_list; params::NamedTuple = (T=291.483, P=1013.2))
    param_names = join(String.(keys(params)), ",")
    ps = "($(param_names))"

    # now replace each of our predefined reaction rates
    rate_list_out = [rate*ps for rate ∈ rate_list]
    rate_replace = Pair.(rate_list, rate_list_out)   # map variable to its initial value
    res = replace(rxnrate_string, rate_replace...)

    res = replace(res,
                  ";" => "",
                  "TEMP" => "T",
                  "EXP" => "exp",
                  "LOG10" => "log10",
                  "@" => "^",
                  "**" => "^",
                  "D-" => "e-",
                  "D+" => "e",
                  "D8" => "e8",
                  "D9" => "e9",
                  "D10" => "e10",
                  "D11" => "e11",
                  "D12" => "e12",
                  "D13" => "e13",
                  "D14" => "e14",
                  "D15" => "e15",
                  "D16" => "e16",
                  "D17" => "e17",
                  "D18" => "e18",
                  "D19" => "e19",
                  "J<" => "J_",
                  "J <" => "J_",  # weird space for some reason
                  ">" => "(t)",
                  # replace number densities by their functions as well
                  "*M" => "*M$(ps)",
                  "*O2" => "*O2$(ps)",
                  "*N2" => "*N2$(ps)",
                  "*H2O" => "*H2O$(ps)",
                  )
    return res
end



# function generate_rrates_mechanism(fac_dict, rate_list; model_name::String="mcm", params::NamedTuple = (T=291.483, P=1013.2))
#     param_names = join(String.(keys(params)), ",")

#     # if file already exists, delete it
#     outpath = "./model/$(model_name)/rrates_mechanism.jl"
#     if isfile(outpath)
#         rm(outpath)
#     end

#     if !isdir("./model/$(model_name)")
#         mkdir("./model/$(model_name)")
#     end

#     species, reactions = parse_rxns(fac_dict["reaction_definitions"])
#     unique_species = [spec for spec ∈ species if spec != nothing]

#     # write the array of reaction to the file
#     open(outpath, "w") do f
#         for i ∈ 1:length(reactions)
#             reactants, reactants_stoich, products, products_stoich, rrate_string = reactions[i]
#             rrate_string = "k_$(i)($(param_names),RO2,t) = "*strip(rxnrate_to_expression(rrate_string, rate_list; params=params))
#             println(f, rrate_string)
#         end
#     end
# end




# NOTE: we can probably do some parallelization here...

# compute_v = """
#     @inbounds for j ∈ axes(R,2)
#         v[j] = k_rates[j]
#         @inbounds for i ∈ nzrange(R, j)
#             v[j] *= u[Rrows[i]]^Rvals[i]
#         end
#     end
# """


compute_v = """

    du .= 0.0

    for i ∈ axes(R,1)
        for k ∈ derivative_terms[i].ks
            du_temp = N[i,k]*k_rates[k]
            for ℓ ∈ reaction_terms[k].is
                du_temp *= u[ℓ]^R[ℓ,k]
            end
            du[i] += du_temp
        end
    end


"""



# update_J = """
#     # update Vmat
#     @inbounds for j ∈ axes(R,2)
#         @inbounds for n ∈ nzrange(R,j)
#             Vmat[j,Rrows[n]] = k_rates[j]
#             @inbounds for i ∈ nzrange(R,j)
#                 if Rrows[i] == Rrows[n]
#                     Vmat[j, Rrows[n]] *= Rvals[i]*u[Rrows[i]]^(Rvals[i]-1)
#                 else
#                     Vmat[j, Rrows[n]] *= u[Rrows[i]]^Rvals[i]
#                 end
#             end
#         end
#     end

#     # update J
#     mul!(Jac, N, Vmat)

# """

update_J = """

    for jac_term ∈ jac_terms
        Jac[jac_term.i,jac_term.j] = 0
        for rxn_term ∈ jac_term.ks
            Jtemp = N[jac_term.i,rxn_term.k] * k_rates[rxn_term.k] * R[jac_term.j,rxn_term.k] * u[jac_term.j]^(R[jac_term.j,rxn_term.k]-1)
            for ℓ ∈ rxn_term.is[2:end]
                Jtemp *= u[ℓ]^R[ℓ,rxn_term.k]
            end
            Jac[jac_term.i,jac_term.j] += Jtemp
        end
    end


"""


function generate_rrates_mechanism(fac_dict, rate_list; model_name::String="mcm", params::NamedTuple = (T=291.483, P=1013.2))
    param_names = join(String.(keys(params)), ",")

    # if file already exists, delete it
    outpath = "./model/$(model_name)/rrates_mechanism.jl"
    if isfile(outpath)
        rm(outpath)
    end

    if !isdir("./model/$(model_name)")
        mkdir("./model/$(model_name)")
    end

    species, reactions = parse_rxns(fac_dict["reaction_definitions"])
    unique_species = [spec for spec ∈ species if spec != nothing]

    open(outpath, "w") do f
        # write the RHS function for the odes
        println(f, "function f!(du, u, p, t)")
        println(f, "\t\t# unpack parameters")
        # println(f, "\t\t$(param_names), N, R, Rrows, Rvals, idx_ro2, RO2, k_rates, v, Vmat = p")
        # println(f, "\t\t$(param_names), N, R, derivative_terms, reaction_terms, idx_ro2, RO2, k_rates, du_temp = p")
        println(f, "\t\t$(param_names), N, R, derivative_terms, reaction_terms, idx_ro2, RO2, k_rates, du_temp, jac_terms, Jtemp = p")
        println(f, "\n\t\t# for peroxy-radical sum")
        println(f, "\t\tRO2 = sum(u[idx_ro2])")

        println(f,"\n\n\t\t# update reaction rate coefficients\n")
        for i ∈ 1:length(reactions)
            reactants, reactants_stoich, products, products_stoich, rrate_string = reactions[i]
            rrate_string = "\t\tk_rates[$(i)]="*strip(rxnrate_to_expression(rrate_string, rate_list; params=params))
            println(f, rrate_string)
        end

        println(f, "\n\t\t# compute reaction rate vector")
        println(f, compute_v)

        #println(f, "\t\t# update RHS with stoich matrix")
        #println(f, "\t\tmul!(du, N, v)")

        println(f, "end")



        # write the Jacobian of the RHS
        println(f, "\n\nfunction Jac!(Jac,u,p,t)")
        println(f,"\t\t# unpack parameters")
        #println(f, "\t\t$(param_names), N, R, Rrows, Rvals, idx_ro2, RO2, k_rates, v, Vmat = p")
        println(f, "\t\t$(param_names), N, R, derivative_terms, reaction_terms, idx_ro2, RO2, k_rates, du_temp, jac_terms, Jtemp = p")
        println(f, "\t\t# for peroxy-radical sum")
        println(f, "\t\tRO2 = sum(u[idx_ro2])")
        println(f,"\n\n\t\t# update reaction rate coefficients\n")
        for i ∈ 1:length(reactions)
            reactants, reactants_stoich, products, products_stoich, rrate_string = reactions[i]
            rrate_string = "\t\tk_rates[$(i)]="*strip(rxnrate_to_expression(rrate_string, rate_list; params=params))
            println(f, rrate_string)
        end

        println(f, update_J)
        println(f, "end")
    end

end


