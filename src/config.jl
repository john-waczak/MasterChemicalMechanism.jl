config_string_old = """

const latitude = 32.949589 # degrees
const longitude = -96.823666  # degrees
const altitude = 0.0
const t_start = DateTime(2023, 3, 7, 6, 0, 0)


# make declarations using parameters for better printing of system

@constants kb = 1.380649e−23 # J/K
@parameters T P

@variables t M O2 N2 H2O  # t is simulation time in minutes

press_pa = 100.0 * P  # 100 Pa / mbar
# 1 Pa = 1 J / m3. We want cm^3 so convert:
press_final = press_pa * 1.0e-6 # there are (10^2)^3 = 10^6 cm³/m³
M = press_final/(kb*T)  # we now have a stand in for number density in molecules/cm³jk

# define O2 and N2 with usual atmospheric abundances
O2 = 0.2095 * M
N2 = 0.7808 * M

# http://www.atmo.arizona.edu/students/courselinks/fall16/atmo336/lectures/sec1/composition.html#:~:text=Water%20vapor%20is%20literally%20individual,0%25%20in%20cold%20polar%20regions
H2O = 0.004 * M # ^ suggest average of 0.4% but between 0-4% globally


# http://www.atmo.arizona.edu/students/courselinks/fall16/atmo336/lectures/sec1/composition.html#:~:text=Water%20vapor%20is%20literally%20individual,0%25%20in%20cold%20polar%20regions
# we should be able to get H2O from water vapor pressure at specific T and P
# e.g. this link:https: //physics.stackexchange.com/questions/596684/what-is-the-formula-to-calculate-the-amount-of-water-in-the-air-at-a-certain-tem


"""


config_string = """


const latitude = 32.949589 # degrees
const longitude = -96.823666  # degrees
const altitude = 0.0
const t_start = DateTime(2023, 3, 7, 6, 0, 0)

const kb = 1.380649e−23 # J/K

function M(T,P)
    press_pa = 100.0 * P  # 100 Pa / mbar
    # 1 Pa = 1 J / m3. We want cm^3 so convert:
    press_final = press_pa * 1.0e-6 # there are (10^2)^3 = 10^6 cm³/m³
    Mout = press_final/(kb*T)  # we now have a stand in for number density in molecules/cm³jk
    return Mout
end


O2(T,P) = 0.2095 * M(T,P)
N2(T,P) = 0.7809 * M(T,P)
H2O(T,P) = 0.004 * M(T,P)

# test it out... These will be our parameters later
params = (
    T = 291.483, # 65 °F
    P = 1013.2  # milibars standard pressure
)

dt = 15.0  # minutes
tspan = (0.0, 60.0)

"""



function generate_config(fac_dict::Dict; model_name::String="mcm")
    # if file already exists, delete it
    outpath = "./model/$(model_name)/config.jl"
    if isfile(outpath)
        rm(outpath)
    end

    if !isdir("./model/$(model_name)")
        mkdir("./model/$(model_name)")
    end

    open(outpath, "a") do f
        println(f, config_string)
    end
end

