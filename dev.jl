# run this file using the project in root of this directory
using ModelingToolkit
using Catalyst
using SolarGeometry
using Dates
using Latexify

using MasterChemicalMechanism
import MasterChemicalMechanism: read_fac_file, parse_rxns

# 1. parse .fac file for text with reaction info
# include("parse_fac.jl")

#fpath = "./src/data/extracted/full/mcm_subset.fac"
#fpath = "./src/data/extracted/monoterpines/limonene.fac"
fpath = "./src/data/extracted/monoterpines/alpha_pinene.fac"
#fpath = "./src/data/extracted/alkanes/methane.fac"
ispath(fpath)

fac_dict = read_fac_file(fpath)
rxns = fac_dict["reaction_definitions"]
species, reactions = parse_rxns(rxns)


# 2. establish simulation parameters, variables
#include("config.jl")

# time, temperature, pressure
@variables t  # simulation time in minutes

@parameters T
@parameters P
# @variables H2O O2 N2


# 3. establish constants
#include("constants.jl")  # e.g. boltzman constant
@constants kb = 1.380649e−23 # J/K

# 4. define relevant atmospheric functions
#include("atmospheric_functions.jl")

press_pa = 100.0 * P  # 100 Pa / mbar
# 1 Pa = 1 J / m3. We want cm^3 so convert:
press_final = press_pa * 1.0e-6 # there are (10^2)^3 = 10^6 cm³/m³
M = press_final/(kb*T)  # we now have a stand in for number density in molecules/cm³jk

# define O2 and N2 with usual atmospheric abundances
O2 = 0.2095 * M
N2 = 0.7808 * M

# http://www.atmo.arizona.edu/students/courselinks/fall16/atmo336/lectures/sec1/composition.html#:~:text=Water%20vapor%20is%20literally%20individual,0%25%20in%20cold%20polar%20regions
H2O = 0.004 * M # ^ suggest average of 0.4% but between 0-4% globally

# we should be able to get H2O from water vapor pressure at specific T and P
# e.g. this link:https: //physics.stackexchange.com/questions/596684/what-is-the-formula-to-calculate-the-amount-of-water-in-the-air-at-a-certain-tem



# 5. generate reaction rate expressions
#include("reaction_rates.jl")

keys(fac_dict)
grcs = fac_dict["generic_rate_coefficients"]
grcs[1]


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

for grc ∈ grcs
    res = grc_to_expression(grc)
    println(res)
    eval(Meta.parse(res))
end


crcs = fac_dict["complex_rate_coefficients"]

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

for crc ∈ crcs
    res = crc_to_expression(crc)
    println(res)
    eval(Meta.parse(res))
end



# 6. define photolysis rates
# include("photolysis.jl")
latitude = 32.949589 # degrees
longitude = -96.823666  # degrees
altitude = 0.0 # meters
t_start = DateTime(2023, 3, 7, 6, 0, 0)

# so this works fine...
res = solar_azimuth_altitude(t_start, latitude, longitude, altitude)


# now we need a function to generate datetime from time in seconds. Closest second should be good enough
# since this is only for box model
test_t = 1.111
t_start + Dates.Second(round(Int, test_t))
Dates.Second(1)

# NOTE: simulation time is in minutes
solar_azimuth(t) = solar_azimuth_altitude(t*60.0, t_start, latitude, longitude, altitude)[1]
solar_altitude(t) = solar_azimuth_altitude(t*60.0, t_start, latitude, longitude, altitude)[2]
# solar zenith angle is complement to solar altitude (aka solar elevation)
SZA(t) = 90 - solar_altitude(t)


function cosx(t)
    res = max(0.0, cosd(SZA(t)))
    return res
end

secx(t) = 1.0/(cosx(t) + eps(typeof(t)))

secx(1.0)

# okay, now we have what we need to define our photolysis rates
# definitions from:  http://chmlin9.leeds.ac.uk/MCMv3.3.1/parameters/photolysis_param.htt

# Inorganics
J_1(t) =6.073E-05*cosx(t)^(1.743)*exp(-1.0*0.474*secx(t))
J_2(t) =4.775E-04*cosx(t)^(0.298)*exp(-1.0*0.080*secx(t))
J_3(t) =1.041E-05*cosx(t)^(0.723)*exp(-1.0*0.279*secx(t))
J_4(t) =1.165E-02*cosx(t)^(0.244)*exp(-1.0*0.267*secx(t))
J_5(t) =2.485E-02*cosx(t)^(0.168)*exp(-1.0*0.108*secx(t))
J_6(t) =1.747E-01*cosx(t)^(0.155)*exp(-1.0*0.125*secx(t))
J_7(t) =2.644E-03*cosx(t)^(0.261)*exp(-1.0*0.288*secx(t))
J_8(t) =9.312E-07*cosx(t)^(1.230)*exp(-1.0*0.307*secx(t))

# Carbonyls
J_11(t) =4.642E-05*cosx(t)^(0.762)*exp(-1.0*0.353*secx(t))
J_12(t) =6.853E-05*cosx(t)^(0.477)*exp(-1.0*0.323*secx(t))
J_13(t) =7.344E-06*cosx(t)^(1.202)*exp(-1.0*0.417*secx(t))
J_14(t) =2.879E-05*cosx(t)^(1.067)*exp(-1.0*0.358*secx(t))
J_15(t) =2.792E-05*cosx(t)^(0.805)*exp(-1.0*0.338*secx(t))
J_16(t) =1.675E-05*cosx(t)^(0.805)*exp(-1.0*0.338*secx(t))
J_17(t) =7.914E-05*cosx(t)^(0.764)*exp(-1.0*0.364*secx(t))
J_18(t) =1.482E-06*cosx(t)^(0.396)*exp(-1.0*0.298*secx(t))
J_19(t) =1.482E-05*cosx(t)^(0.396)*exp(-1.0*0.298*secx(t))
J_20(t) =7.600E-04*cosx(t)^(0.396)*exp(-1.0*0.298*secx(t))
J_21(t) =7.992E-07*cosx(t)^(1.578)*exp(-1.0*0.271*secx(t))
J_22(t) =5.804E-06*cosx(t)^(1.092)*exp(-1.0*0.377*secx(t))
J_23(t) =2.4246E-06*cosx(t)^(0.395)*exp(-1.0*0.296*secx(t))
J_24(t) =2.424E-06*cosx(t)^(0.395)*exp(-1.0*0.296*secx(t))
J_31(t) =6.845E-05*cosx(t)^(0.130)*exp(-1.0*0.201*secx(t))
J_32(t) =1.032E-05*cosx(t)^(0.130)*exp(-1.0*0.201*secx(t))
J_33(t) =3.802E-05*cosx(t)^(0.644)*exp(-1.0*0.312*secx(t))
J_34(t) =1.537E-04*cosx(t)^(0.170)*exp(-1.0*0.208*secx(t))
J_35(t) =3.326E-04*cosx(t)^(0.148)*exp(-1.0*0.215*secx(t))

# Organic Peroxides
J_41(t) =7.649E-06*cosx(t)^(0.682)*exp(-1.0*0.279*secx(t))

# Organic Nitrates
J_51(t) =1.588E-06*cosx(t)^(1.154)*exp(-1.0*0.318*secx(t))
J_52(t) =1.907E-06*cosx(t)^(1.244)*exp(-1.0*0.335*secx(t))
J_53(t) =2.485E-06*cosx(t)^(1.196)*exp(-1.0*0.328*secx(t))
J_54(t) =4.095E-06*cosx(t)^(1.111)*exp(-1.0*0.316*secx(t))
J_55(t) =1.135E-05*cosx(t)^(0.974)*exp(-1.0*0.309*secx(t))
J_56(t) =4.365E-05*cosx(t)^(1.089)*exp(-1.0*0.323*secx(t))
J_57(t) =3.363E-06*cosx(t)^(1.296)*exp(-1.0*0.322*secx(t))
J_61(t) =7.537E-04*cosx(t)^(0.499)*exp(-1.0*0.266*secx(t))


times = 0.0:1.0:(60*60*24)
j1 = J_1.(times)


# 7. Generate list of species

unique_species = [spec for spec ∈ species if spec != nothing]
length(unique_species)
@species (X(t))[1:length(unique_species)]  # these will be our concentrations


# to replace, we simply look up the idx of the species and use the corresponding array element
# e.g.
# unique_species[21]
# idx = findfirst(elem -> elem == "CH3OH", unique_species)
# X[21]


# 8. Peroxy Radical Sum
ro2_string = fac_dict["peroxy_radicals"]

function get_RO2_species(ro2_string, unique_species)
    ro2_string = replace(ro2_string, ";"=>"")
    ro2_sum = split(ro2_string, "=")[2]
    ro2_species = [strip(spec) for spec ∈ split(ro2_sum, "+")]
    ro2_symbolic = []
    for spec ∈ ro2_species
        idx = findfirst(elem -> elem == spec, unique_species)
        push!(ro2_symbolic, X[idx])
    end

    return ro2_symbolic
end

ro2_species = get_RO2_species(ro2_string, unique_species)

length(ro2_species)

RO2 = sum(ro2_species)
RO2


# 9. evaluate reaction rate coefficients
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

test_rxn_rate = reactions[1][end]
rxnrate_to_expression(test_rxn_rate)

# 10. Generate reactions


# need a function for identifying repeats
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


test_species = ["H2O", "H2O", "N2O"]
test_stoich = [2.0, 1.0, 1.0]

out_species, out_stoich = remove_duplicates(test_species, test_stoich)

@assert out_species == unique(test_species)
@assert all(out_stoich .== [3.0, 1.0])


rxns = []
for rxn ∈ reactions
    reactants, reactants_stoich, products, products_stoich, rrate_string = rxn

    # make sure stoichiometric coefficients are ints
    reactants_stoich = [Int(rs) for rs ∈ reactants_stoich]
    products_stoich = [Int(ps) for ps ∈ products_stoich]

    # check for duplicates and update stoich
    reactants, reactants_stoich = remove_duplicates(reactants, reactants_stoich)
    products, products_stoich = remove_duplicates(products, products_stoich)

    eval(Meta.parse("rrate = "*rxnrate_to_expression(rrate_string)))

    in = []
    if reactants == [nothing]
        in = nothing
        reactants_stoich = nothing
    else
        for reactant ∈ reactants
            idx = findfirst(x -> x == reactant, unique_species)
            push!(in, X[idx])
        end
    end

    out = []
    if products == [nothing]
        out = nothing
        products_stoich = nothing
    else
        for product ∈ products
            idx = findfirst(x -> x == product, unique_species)
            push!(out, X[idx])
        end
    end

    push!(rxns, Reaction(rrate,
                         in,
                         out,
                         reactants_stoich,
                         products_stoich
                         ))
end


reactions[1]
rxns[1]

# set simulation timescale to minutes
@named mcm = ReactionSystem(rxns, t, collect(X), [T, P])

println(parameters(mcm))
# println(reactions(mcm))


# too much memory for tex to handle
# ltx = latexify(mcm)
# # save to file
# open("mechanism.qmd", "w") do f
#     println(f, ltx)
# end


using DifferentialEquations

ndays = 25
tspan = (0.0, ndays*24.0*60.0)

# use these initial concentrations and set everything
# else to 0.0: https://github.com/AtChem/AtChem2/blob/master/model/configuration/initialConcentrations.config



# init_dict = Dict(
#     "CH4" => 4.9e+13,
#     "CO" => 3.6e+12,
#     "O3" => 5.2e+11,
#     "NO2" => 2.4e+11,
# )

m_init = substitute(M, Dict([
    T => 291.483, # 65 °F
    P => 1013.2,  # milibars standard pressure
    kb => 1.380649e−23 # J/K
]))
m_init = m_init.val
nₘ = 10^9 # for ppb

# convert ppb to molecules/cc
init_dict = Dict(
    "CH4" => 1909.6*m_init/nₘ,
    "CO" => 100.0*m_init/nₘ,
    "O3" => 18.0*m_init/nₘ,
    "NO2" => 100*m_init/nₘ,
    "APINENE" => 30.0*m_init/nₘ,
    "H2" => 500*m_init/nₘ,
)


# see this link: https://docs.sciml.ai/Catalyst/stable/example_networks/smoluchowski_coagulation_equation/#smoluchowski_coagulation_equation
u₀    = zeros(Float64, size(X,1))
for (key, val) ∈ init_dict
    try
        println("$(key): $(val)")
        idx = findfirst(x -> x == key, unique_species)
        println("idx: ", idx)
        u₀[idx] = val
    catch e
        println(e)
    end
end
u₀map = Pair.(collect(X), u₀)   # map variable to its initial value

params = (
    T => 291.483, # 65 °F
    P => 1013.2  # milibars standard pressure
)




# let's convert the model into an ODESystem so we can simplify

odesys = convert(ODESystem, mcm; combinatoric_ratelaws=false)
#ode_simplified = structural_simplify(odesys)
#ode_prob = ODEProblem(ode_simplified, u₀, tspan, params; jac=true, sparse=true)
# ode_prob = ODEProblem(odesys, u₀, tspan, params; jac=true, sparse=true)
ode_prob = ODEProblem(odesys, u₀, tspan, params; sparse=true)


sol = solve(ode_prob, saveat=15)


# got it to work for the limonene system
println(numspecies(mcm))  # 3490
println(numreactions(mcm)) # 1163

println(unique_species)

using Plots

unique_species

speciesiwant = ["CO", "O3", "NO2", "APINENE", "H2"]
idxs= [findfirst(x->x==spec, unique_species) for spec ∈ speciesiwant]

p = plot()
i = 1
for idx ∈ idxs
    plot!(p, sol.t ./ (60*24), sol[idx, :] .* (nₘ/m_init), label=speciesiwant[i], xlabel="t [days]", ylabel="concentration [ppb]", legend=:outertopright)
    i += 1
end


display(p)


