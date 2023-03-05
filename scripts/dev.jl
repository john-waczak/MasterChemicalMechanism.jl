
fpath = "../src/data/extracted/full/mcm_subset.fac"
isfile(fpath)

fac_file = String(read(fpath))

function read_fac_file(path::String)
    out_dict = Dict()

    fac_file = String(path)
    split_file = split(fac_file, "*;")

    out_dict["VOCs"] = split(split_file[3], "\r\n")[3:end-2]
    out_dict["VOCs"] = [elem for elem ∈ split(join(out_dict["VOCs"], " "), " ") if !(elem ∈ (";", ""))]

    out_dict["generic_rate_coefficients"] = split(split_file[5], "\r\n")[2:end-1]

    out_dict["complex_rate_coefficients"] = split(split_file[7], "\r\n")[2:end-2]


    # turn this into a julia expression we can write to a file
    out_dict["peroxy_radicals"] = join(lstrip.(split(split_file[9], "\r\n")[6:end-1]), " ")

    # out_dict["reaction_definitions"] = split(split_file[11], "\r\n")[2:end-1]
    # use regex instead
    eqn_idxs = findall(r"%(.*?)\;",fac_file)  # find anything between } and ;
    out_dict["reaction_definitions"] = [fac_file[eqn_idx] for eqn_idx ∈ eqn_idxs]


    return out_dict
end

fac_dict = read_fac_file(fac_file)

fac_dict["VOCs"]
fac_dict["generic_rate_coefficients"];
fac_dict["generic_rate_coefficients"][1]
fac_dict["generic_rate_coefficients"][end]

fac_dict["complex_rate_coefficients"];
fac_dict["complex_rate_coefficients"][1]
fac_dict["complex_rate_coefficients"][end]

fac_dict["peroxy_radicals"]

fac_dict["reaction_definitions"]



rxns = fac_dict["reaction_definitions"]
rxns = [strip(replace(rxn, "%"=>"", ";"=>"", "\""=>"")) for rxn ∈ rxns]


function species_and_stoich(rxn_half)
    species_list = strip.(split(rxn_half, "+"))

    involved_species = []
    species_stoich = Float32[]

    if species_list == [""]
        push!(involved_species, nothing)
        push!(species_stoich, 1.0)
    else
        for species ∈ species_list
            # find index of first non-number (i.e. ignore stoich)
            idx = findfirst(isletter.(collect(species)))
            if idx == 1
                stoich = 1.0
                species = species
            else
                stoich = parse(Float32, species[1:idx-1])
                species = species[idx:end]
            end

            push!(involved_species, species)
            push!(species_stoich, stoich)
        end
    end

    return involved_species, species_stoich
end


function parse_rxns(rxns)
    species_list = []
    reactions = []

    for rxn ∈ rxns
        # separate rate and reaction
        rate_equation, reaction_full = split(rxn, ":")

        # split into reactants and products
        lhs, rhs = split(reaction_full, "=")

        # parse reaction info
        reactants, reactants_stoich = species_and_stoich(lhs)
        products, products_stoich = species_and_stoich(rhs)

        push!(species_list, reactants)
        push!(species_list, products)
        push!(reactions, (reactants, reactants_stoich, products, products_stoich, rate_equation))
    end

    return unique(vcat(species_list...)), reactions
end

species, reactions = parse_rxns(rxns)

size(species)
size(reactions)


# m=2.55e+19
# o2=0.2095*m


# REAL(dp)::M, N2, O2, RO2, H2O, TEMP,

# We should specify
# temp, pressure, H2O as variables which we can set initially or read from datafile


# RO2 list: https://github.com/AtChem/AtChem2/blob/master/mcm/peroxy-radicals_v3.3.1
# photolysis rates: https://github.com/AtChem/AtChem2/blob/master/mcm/photolysis-rates_v3.3.1

# atmospheric functions (i.e. M calculation) https://github.com/AtChem/AtChem2/blob/37503aefb02bb54fa3b04899c68d15afa8b1c8c0/src/atmosphereFunctions.f90
