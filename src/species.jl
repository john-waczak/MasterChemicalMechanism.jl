# function generate_species(fac_dict::Dict model_name::String="mcm")
#     # if file already exists, delete it
#     outpath = "./model/$(model_name)/species.jl"
#     if isfile(outpath)
#         rm(outpath)
#     end

#     if !isdir("./model/$(model_name)")
#         mkdir("./model/$(model_name)")
#     end



#     rxns = fac_dict["reaction_definitions"]
#     species, reactions = parse_rxns(rxns)
#     unique_species = [spec for spec ∈ species if spec != nothing]

#     open(outpath, "w") do f
#         println(f, "# -------------------------")
#         println(f, "# species list")
#         println(f, "# -------------------------")
#         println(f,"species = [")
#         for spec ∈ unique_species
#             println(f, "\t\"$(spec)\",")
#         end
#         println(f,"]\n\n")

#         println(f, "@species (X(t))[1:length(species)]")
#         println("There are $(length(unique_species)) unique species in this mechanism")
#     end

# end


function generate_species(fac_dict::Dict; model_name::String="mcm")
    # if file already exists, delete it
    outpath = "./model/$(model_name)/species.jl"
    if isfile(outpath)
        rm(outpath)
    end

    if !isdir("./model/$(model_name)")
        mkdir("./model/$(model_name)")
    end



    rxns = fac_dict["reaction_definitions"]
    species, reactions = parse_rxns(rxns)
    unique_species = [spec for spec ∈ species if spec != nothing]

    open(outpath, "w") do f
        println(f, "# -------------------------")
        println(f, "# species list")
        println(f, "# -------------------------")
        println(f,"species = [")
        for spec ∈ unique_species
            println(f, "\t\"$(spec)\",")
        end
        println(f,"]\n\n")
        println("There are $(length(unique_species)) unique species in this mechanism")
    end

end
