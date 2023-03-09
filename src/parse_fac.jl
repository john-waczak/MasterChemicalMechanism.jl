"""
    function read_fac_file(path::String)

Parse a mechanism file in the FACSIMILE format (ending with .fac). This assumes that the file has all necessary sections to define the master chemical mechanism including

- `VOCs`: a list of VOCs tracked in the simulation
- `generic_rate_coefficients`
- `complex_rate_coefficients`
- `peroxy_radicals`: the species being summed to form `RO2`
- `reaction_definitions`: The list of chemical reactions in our mechanism

# Output
Returns a dictionary containing the above sections as text.
"""
function read_fac_file(path::String)
    out_dict = Dict()

    fac_file = String(read(path))
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
    out_dict["reaction_definitions"] = [strip(replace(rxn, "%"=>"", ";"=>"", "\""=>"")) for rxn ∈ out_dict["reaction_definitions"]]

    return out_dict
end



"""
    function species_and_stoich(rxn_half)

Given half of a reaction (either reactants or products), parse the equation and return a list of species involved and associated stoichiometric coefficients. If no species is involved (production or destruction), set the species to `nothing`.
"""
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



"""
    function parse_rxns(rxns)

Given a list of equations `rxns`, parse equation returning a list of unique species in the entire mechanism and a list of reaction components, i.e. tuples of `(reactants, reactants_stoich, products, products_stoich, rate_equation)`
"""
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
