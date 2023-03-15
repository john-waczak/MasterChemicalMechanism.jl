function grc_to_expression(grc_string)
    return replace(grc_string,
                   ";" => "",
                   "TEMP" => "T",
                   "EXP" => "exp",
                   "D-" => "e-",
                   "D+" => "e",
                   "D7" => "e7"
                   )
end

function crc_to_expression(crc_string)
    return replace(crc_string,
                   ";" => "",
                   "TEMP" => "T",
                   "EXP" => "exp",
                   "LOG10" => "log10",
                   "@" => "^",
                   "**" => "^",
                   "D-" => "e-",
                   "D+" => "e",
                   "D17" => "e17",
                   )
end


function generate_rrates(fac_dict::Dict, model_name::String="mcm")
    # if file already exists, delete it
    outpath = "./model/$(model_name)/rrates.jl"
    if isfile(outpath)
        rm(outpath)
    end

    if !isdir("./model/$(model_name)")
        mkdir("./model/$(model_name)")
    end


    # write generic rate coefficients
    grcs = fac_dict["generic_rate_coefficients"]
    open(outpath, "a") do f
        println(f, "# -------------------------")
        println(f, "# generic rate coefficients")
        println(f, "# -------------------------")

        for grc ∈ grcs
            res = grc_to_expression(grc)
            println(res)
            println(f, rstrip(res))
        end

        println(f, "\n\n")
    end


    # write complex rate coefficients
    crcs = fac_dict["complex_rate_coefficients"]
    open(outpath, "a") do f
        println(f, "# -------------------------")
        println(f, "# generic rate coefficients")
        println(f, "# -------------------------")
        for crc ∈ crcs
            res = crc_to_expression(crc)
            println(res)
            println(f, rstrip(res))
        end
    end
end

