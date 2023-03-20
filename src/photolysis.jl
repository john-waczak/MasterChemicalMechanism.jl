mcm_photolysis_model = """
using SolarGeometry

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
"""


function generate_photolysis_mcm(fac_dict::Dict; model_name::String="mcm")
    # if file already exists, delete it
    outpath = "./model/$(model_name)/photolysis.jl"
    if isfile(outpath)
        rm(outpath)
    end

    if !isdir("./model/$(model_name)")
        mkdir("./model/$(model_name)")
    end

    open(outpath, "a") do f
        println(f, mcm_photolysis_model)
    end
end
