# structs for configuring simulation

struct SimulationConfig
    mechanism_file::String
    temp::typeof(1.0*u"K")
    RH::Float64  # relative humidity ∈ [0, 1.0]
    start_time::DateTime
    end_time::DateTime
    H2O::typeof(1.0*u"1/cm^3") # water concentration
    latitude::typeof(1.0*u"°")
    longitude::typeof(1.0*u"°")
end
