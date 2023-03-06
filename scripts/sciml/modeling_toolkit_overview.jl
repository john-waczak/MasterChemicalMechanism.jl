using ModelingToolkit
using Plots, DifferentialEquations


@variables t x(t)    # independent and dependent variables
@parameters τ        # parameters
@constants h = 1     # constants have an assigned value

D = Differential(t)  # define an operator for differentiation w.r.t. time

# define an ODE via symbolic equation using  ~ (tilde) for equality since = is assignment operator
@named fol_model = ODESystem([
    D(x) ~ (h-x) / τ
])

prob = ODEProblem(fol_model, [x => 0.0], (0.0, 10.0), [τ => 3.0])  # parameter τ can be given a value but h can't

sol = solve(prob)

plot(sol)



# algebraic relations and structural simplification
@variables RHS(t)
@named fol_separate = ODESystem([
    RHS ~ (h-x) / τ,
    D(x) ~ RHS
])

fol_simplified = structural_simplify(fol_separate)
equations(fol_simplified)

# ok so our simplified model is equivalent to the original model
equations(fol_simplified) == equations(fol_model)
