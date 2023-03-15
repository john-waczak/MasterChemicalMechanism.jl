
function rxnrate_to_expression(rxnrate_string)
      return replace(rxnrate_string,
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
                   )
end



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




function generate_mechanism(fac_dict::Dict, model_name::String="mcm")
    # if file already exists, delete it
    outpath = "./model/$(model_name)/mechanism.jl"
    if isfile(outpath)
        rm(outpath)
    end

    if !isdir("./model/$(model_name)")
        mkdir("./model/$(model_name)")
    end

    species, reactions = parse_rxns(fac_dict["reaction_definitions"])
    unique_species = [spec for spec ∈ species if spec != nothing]

    # write the array of reaction to the file
    open(outpath, "w") do f
        println(f, "rxns = [")
        for rxn ∈ reactions
            reactants, reactants_stoich, products, products_stoich, rrate_string = rxn

            # make sure stoichiometric coefficients are ints
            reactants_stoich = [Int(rs) for rs ∈ reactants_stoich]
            products_stoich = [Int(ps) for ps ∈ products_stoich]

            # check for duplicates and update stoich
            reactants, reactants_stoich = remove_duplicates(reactants, reactants_stoich)
            products, products_stoich = remove_duplicates(products, products_stoich)

            rrate_string = rxnrate_to_expression(rrate_string)
            # eval(Meta.parse("rrate = "*rxnrate_to_expression(rrate_string)))

            in = "["
            if reactants == [nothing]
                in = nothing
                reactants_stoich = nothing
            else
                for reactant ∈ reactants
                    idx = findfirst(x -> x == reactant, unique_species)
                    in = in*"X[$(idx)], "
                end
                in = in * "]"
            end

            out = "["
            if products == [nothing]
                out = nothing
                products_stoich = nothing
            else
                for product ∈ products
                    idx = findfirst(x -> x == product, unique_species)
                    out = out * "X[$(idx)], "
                end
                out = out * "]"
            end

            println(f, "\t\tReaction($(rrate_string), $(in), $(out), $(reactants_stoich), $(products_stoich)),")

            # push!(rxns, Reaction(rrate,
            #                     in,
            #                     out,
            #                     reactants_stoich,
            #                     products_stoich
            #                     ))
        end

        println(f, "]\n\n")

        # write code to generate named reaction system
        println(f, "@named mechanism = ReactionSystem(rxns, t, collect(X), [T, P])")
    end
end
