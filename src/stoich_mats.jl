using SparseArrays
using CSV, DataFrames

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



function generate_stoich_mat(fac_dict::Dict; model_name::String="mcm")
    # if file already exists, delete it
    outpath = "./model/$(model_name)/N.csv"
    if isfile(outpath)
        rm(outpath)
    end

    outpath2 = "./model/$(model_name)/R.csv"

    if !isdir("./model/$(model_name)")
        mkdir("./model/$(model_name)")
    end

    species, reactions = parse_rxns(fac_dict["reaction_definitions"])
    unique_species = [spec for spec ∈ species if spec != nothing]


    N = zeros(Int, length(unique_species), length(reactions))

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

        # update N
        if idxs_reactants != nothing
            for j ∈ 1:size(idxs_reactants,1)
                N[idxs_reactants[j], i] = -1*reactants_stoich[j]
            end
        end

        if idxs_products != nothing
            for j ∈ 1:size(idxs_products,1)
                N[idxs_products[j], i] = products_stoich[j]
            end
        end

    end

    # sparsify
    N = sparse(N)
    I,J,V = findnz(N)

    # write sparse matrix to csv file for easy loading
    df = DataFrame([:I => I, :J => J, :V => V])
    CSV.write(outpath, df)


    # now we want to generate our reactants stoich matrix
    R = copy(N)
    R[R .> 0 ] .= 0
    R = -1 .* R
    I,J,V = findnz(R)

    df = DataFrame([:I => I, :J => J, :V => V])
    CSV.write(outpath2, df)
end



function get_sparse_mat(;mat_name::String="N", model_name::String="mcm")
    fpath = "./model/$(model_name)/$(mat_name).csv"
    @assert isfile(fpath)

    df = CSV.File(fpath) |> DataFrame
    Mat = sparse(df.I, df.J, df.V)

    return Mat
end



# function get_stoich_mats(;model_name::String="mcm")
#     fpath = "./model/$(model_name)/N.csv"
#     @assert isfile(fpath)

#     df = CSV.File(fpath) |> DataFrame

#     N = sparse(df.I, df.J, df.V)

#     # now we want to generate our reactants stoich matrix
#     R = copy(N)
#     R[R .> 0 ] .= 0
#     R = -1 .* R

#     return N, R
# end
