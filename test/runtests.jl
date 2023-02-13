using MasterChemicalMechanism
using Test

@testset "atmospheric_functions.jl" begin
    M_stp = MasterChemicalMechanism.calcAirDensity(1013.2, 273.15)
    @test isapprox(M_stp, 2.6866e19, rtol=1e-4)
end
