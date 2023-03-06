# see here: https://github.com/AtChem/AtChem2/blob/37503aefb02bb54fa3b04899c68d15afa8b1c8c0/src/atmosphereFunctions.f90


"""
    calcAirDensity(press, temp)

Calculate the number density of air in [molecules cm-3] from the pressure in [mbar] and the temperature in [K].
"""
function calcAirDensity(press, temp)

    # ----- math -----------------
    # PV = N*kb*T
    # --> n = N/V = P/(kb*T)
    # ----------------------------

    press_pa = press * 100.0  # 100 Pa / mbar

    # 1 Pa = 1 J / m3. We want cm^3 so convert:
    press_final = press_pa * 1.0e-6  # there are (10^2)^3 = 10^6 cm³/m³

    kb = 1.380649e−23 # J/K

    return press_final/(kb*temp)
end


# we should do this symbolically
