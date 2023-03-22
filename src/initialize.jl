function generate_init_dict(init_path, M)
    chem_ignore = ["O2", "N2"]

    init_raw = readlines(init_path)

    out_dict = Dict()
    for line ∈ init_raw
        name, valstring = split(line, "=")
        valstring = replace(valstring,
                            ";"=>"",
                            "D"=>"e",
                            "M"=>"$(M)"
                            )
        valstring = strip(valstring)
        out_dict[name] = eval(Meta.parse(valstring))
    end

    return out_dict
end
