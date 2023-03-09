using MasterChemicalMechanism
using Test

fpath = "../src/data/extracted/full/mcm_subset.fac"
@testset "parse_fac.jl" begin
    import MasterChemicalMechanism: read_fac_file, species_and_stoich, parse_rxns

    @test isfile(fpath)

    fac_dict = read_fac_file(fpath)
    @test Set(keys(fac_dict)) == Set(["VOCs", "generic_rate_coefficients", "complex_rate_coefficients", "peroxy_radicals", "reaction_definitions"])

    rxns = fac_dict["reaction_definitions"]

    test_half = split(split(rxns[1], ":")[2], "=")[1]
    specs, stoichs = species_and_stoich(test_half)
    @test specs[1] == "O"
    @test stoichs[1] == 1.0

    species, reactions = parse_rxns(rxns)
    @test size(species, 1) == 5833
    @test size(reactions, 1) == 17199
end

# see this link for sample mechanism to run tests against:
# https://github.com/AtChem/AtChem2/blob/master/tests/model_tests/env_model_2/env_model_2.fac

