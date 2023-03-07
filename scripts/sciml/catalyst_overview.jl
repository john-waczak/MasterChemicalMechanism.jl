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


# see this link for generating reactions programatically with symbolic array
# https://docs.sciml.ai/Catalyst/dev/example_networks/smoluchowski_coagulation_equation/

# state variables are X, pars stores rate parameters for each rx

@variables t
@parameters α β
@species (X(t))[1:3]

rx = [
    Reaction(α, [X[1], X[2]], [X[2]], [1, 1], [2]),
    Reaction(β, [X[2]], [X[3]], [1], [1])
]

@named rs = ReactionSystem(rx, t, collect(X), [α, β])

p = [:α => 0.1/1000, :β => 0.01]
tspan = (0.0, 250.0)
u0 = [X[1] => 999.0, X[2] => 1.0, X[3] => 0.0]

op = ODEProblem(rs, u0, tspan, p)
sol = solve(op)

plot(sol, lw=2)

