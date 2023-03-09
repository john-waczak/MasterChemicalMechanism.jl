# assuming pressure P is defined,
# see here: https://github.com/AtChem/AtChem2/blob/37503aefb02bb54fa3b04899c68d15afa8b1c8c0/src/atmosphereFunctions.f90


# ----- math -----------------
# PV = N*kb*T
# --> n = N/V = P/(kb*T)
# ----------------------------

press_pa = 100.0 * P  # 100 Pa / mbar
# 1 Pa = 1 J / m3. We want cm^3 so convert:
press_final = press_pa * 1.0e-6 # there are (10^2)^3 = 10^6 cm³/m³

M = press_final/(kb*T)  # we now have a stand in for number density in molecules/cm³

