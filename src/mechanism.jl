
# function rxnrate_to_expression(rxnrate_string)
#       return replace(rxnrate_string,
#                      ";" => "",
#                      "TEMP" => "T",
#                      "EXP" => "exp",
#                      "LOG10" => "log10",
#                      "@" => "^",
#                      "**" => "^",
#                      "D-" => "e-",
#                      "D+" => "e",
#                      "D8" => "e8",
#                      "D9" => "e9",
#                      "D10" => "e10",
#                      "D11" => "e11",
#                      "D12" => "e12",
#                      "D13" => "e13",
#                      "D14" => "e14",
#                      "D15" => "e15",
#                      "D16" => "e16",
#                      "D17" => "e17",
#                      "D18" => "e18",
#                      "D19" => "e19",
#                      "J<" => "J_",
#                      "J <" => "J_",  # weird space for some reason
#                      ">" => "(t)",
#                    )
# end



# function remove_duplicates(species, stoich)
#     unique_species = unique(species)

#     out_species = []
#     out_stoich = []
#     for uspec ∈ unique_species
#         idxs = findall(x -> x == uspec, species)
#         push!(out_species, uspec)
#         push!(out_stoich, sum(stoich[idxs]))
#     end
#     return out_species, out_stoich
# end




# function generate_mechanism(fac_dict::Dict, model_name::String="mcm")
#     # if file already exists, delete it
#     outpath = "./model/$(model_name)/mechanism.jl"
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
#         println(f, "rxns = [")
#         for rxn ∈ reactions
#             reactants, reactants_stoich, products, products_stoich, rrate_string = rxn

#             # make sure stoichiometric coefficients are ints
#             reactants_stoich = [Int(rs) for rs ∈ reactants_stoich]
#             products_stoich = [Int(ps) for ps ∈ products_stoich]

#             # check for duplicates and update stoich
#             reactants, reactants_stoich = remove_duplicates(reactants, reactants_stoich)
#             products, products_stoich = remove_duplicates(products, products_stoich)

#             rrate_string = rxnrate_to_expression(rrate_string)
#             # eval(Meta.parse("rrate = "*rxnrate_to_expression(rrate_string)))

#             in = "["
#             if reactants == [nothing]
#                 in = nothing
#                 reactants_stoich = nothing
#             else
#                 for reactant ∈ reactants
#                     idx = findfirst(x -> x == reactant, unique_species)
#                     in = in*"X[$(idx)], "
#                 end
#                 in = in * "]"
#             end

#             out = "["
#             if products == [nothing]
#                 out = nothing
#                 products_stoich = nothing
#             else
#                 for product ∈ products
#                     idx = findfirst(x -> x == product, unique_species)
#                     out = out * "X[$(idx)], "
#                 end
#                 out = out * "]"
#             end

#             println(f, "\t\tReaction($(rrate_string), $(in), $(out), $(reactants_stoich), $(products_stoich)),")

#             # push!(rxns, Reaction(rrate,
#             #                     in,
#             #                     out,
#             #                     reactants_stoich,
#             #                     products_stoich
#             #                     ))
#         end

#         println(f, "]\n\n")

#         # write code to generate named reaction system
#         println(f, "@named mechanism = ReactionSystem(rxns, t, collect(X), [T, P])")
#     end
# end





function remove_duplicates(species, stoich)
    unique_species = unique(species)

    out_species = []
    out_stoich = []
    for uspec ∈ unique_species
        idxs = findall(x -> x == uspec, species)
        push!(out_species, uspec)
        push!(out_stoich, sum(stoich[idxs]))
    end
    return out_species, out_stoich
end


function get_species_idxs(species, unique_species)
    idxs = []
    for spec ∈ species
        idx = findfirst(x -> x == spec, unique_species)
        push!(idxs, idx)
    end

    if idxs == [nothing]
        idxs = nothing
    end

    return idxs
end


function generate_ode_f(fac_dict::Dict, idx_ro2; model_name::String="mcm", params::NamedTuple = (T=291.483, P=1013.2))
    param_names = join(String.(keys(params)), ",")
    ps = "($(param_names))"

    # if file already exists, delete it
    outpath = "./model/$(model_name)/ode_f.jl"
    if isfile(outpath)
        rm(outpath)
    end

    if !isdir("./model/$(model_name)")
        mkdir("./model/$(model_name)")
    end

    species, reactions = parse_rxns(fac_dict["reaction_definitions"])
    unique_species = [spec for spec ∈ species if spec != nothing]

    du = ["du[$(i)] = " for i ∈ 1:length(unique_species)]

    # loop through reactions and add terms to rhs
    for i ∈ 1:length(reactions)

        reactants, reactants_stoich, products, products_stoich, rrate_string = reactions[i]

        # make sure stoichiometric coefficients are ints
        reactants_stoich = [Int(rs) for rs ∈ reactants_stoich]
        products_stoich = [Int(ps) for ps ∈ products_stoich]

        # check for duplicates and update stoich
        reactants, reactants_stoich = remove_duplicates(reactants, reactants_stoich)
        products, products_stoich = remove_duplicates(products, products_stoich)

        # generate index lookups for reactants and products
        idxs_reactants = get_species_idxs(reactants, unique_species)
        idxs_products = get_species_idxs(products, unique_species)

        # -------------------------------
        # a reaction of the form
        # aA + bB → cC + dD
        # contributes the following terms to the rates
        # d[A]/dt += - a*k*[A]ᵃ*[B]ᵇ
        # d[B]/dt += - b*k*[A]ᵃ*[B]ᵇ
        # d[C]/dt += c*k*[A]ᵃ*[B]ᵇ
        # d[D]/dt += d*k*[A]ᵃ*[B]ᵇ
        # -------------------------------

        # compute reaction rate
        rᵢ = "k_$(i)($(param_names),RO2,t)"
        for j ∈ 1:length(idxs_reactants)
            rᵢ *= "*u[$(idxs_reactants[j])]^$(reactants_stoich[j])"
        end

        # update reactants
        for j ∈ 1:length(idxs_reactants)
            du[idxs_reactants[j]] *= " - $(reactants_stoich[j])*$(rᵢ)"
        end

        # update products
        if idxs_products != nothing
            for j ∈ 1:length(idxs_products)
                du[idxs_products[j]] *= " + $(products_stoich[j])*$(rᵢ)"
            end
        end
    end

    # write the array of reaction to the file
    open(outpath, "w") do f
        println(f, "function mcm!(du, u, p, t)")
        println(f, "\t$(param_names),idx_ro2 = p")
        println(f, "\tRO2 = sum(u[idx_ro2])")
        for eqn ∈ du
            println(f, "\t$(eqn)")
        end
        println(f, "end")
    end
end


