using ModelingToolkit
using DifferentialEquations
using Plots
using Catalyst
using Latexify

rn = @reaction_network begin
    α, S + I --> 2I
    β, I --> R
end

p = [:α => 0.1/1000, :β => 0.01]
tspan = (0.0, 250.0)
u0 = [:S => 999.0, :I => 1.0, :R => 0.0]

op = ODEProblem(rn, u0, tspan, p)
sol = solve(op)

plot(sol, lw=2)
