function generate_ro2(fac_dict::Dict, model_name::String="mcm")
    outpath = "./model/$(model_name)/ro2.jl"

    # if it already exists, remove it so we can recreate it
    if isfile(outpath)
        rm(outpath)
    end

    if !isdir("./model/$(model_name)")
        mkdir("./model/$(model_name)")
    end


    ro2_string = fac_dict["peroxy_radicals"]
    rxns = fac_dict["reaction_definitions"]
    species, reactions = parse_rxns(rxns)
    unique_species = [spec for spec ∈ species if spec != nothing]


    ro2_string = replace(ro2_string, ";"=>"")
    ro2_sum = split(ro2_string, "=")[2]
    ro2_species = [strip(spec) for spec ∈ split(ro2_sum, "+")]

    ro2_symbolic = []
    for spec ∈ ro2_species
        idx = findfirst(elem -> elem == spec, unique_species)
        push!(ro2_symbolic, "X[$(idx)]")
    end

    RO2_string = "RO2 = " * join(ro2_symbolic, "+")
    open(outpath, "w") do f
        println(f, RO2_string)
    end

    println(RO2_string)

    # make another file with the actual species names
    outpath2 = "./model/$(model_name)/ro2_original.jl"
    if isfile(outpath2)
        rm(outpath2)
    end
    open(outpath2, "w") do f
        println(f, ro2_string)
    end

    println(ro2_string)
end

