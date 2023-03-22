using DifferentialEquations
using Plots
using DataFrames, CSV

# standard example

function lorenz!(du, u, p, t)
    du[1] = 10.0*(u[2]-u[1])
    du[2] = u[1]*(28.0-u[3]) - u[2]
    du[3] = u[1]*u[2] - (8/3)*u[3]
end


u0 = [1.0, 0.0, 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz!, u0, tspan)

sol = solve(prob)


plot(sol, idxs=(1,2,3))



# using parametrized functions

function parametrized_lorenz!(du, u, p, t)
    du[1] = p[1]*(u[2]-u[1])
    du[2] = u[1]*(p[2]-u[3]) - u[2]
    du[3] = u[1]*u[2] - p[3]*u[3]
end

u0 = [1.0, 0.0, 0.0]
tspan = (0.0, 100.0)
p = [10.0, 28.0, 8/3]

prob = ODEProblem(parametrized_lorenz!, u0, tspan, p)

sol = solve(prob)
plot(sol, idxs=(1,2,3))

# alternatively
function parametrized_lorenz!(du, u, p, t)
    x,y,z = u
    σ,ρ,β = p
    du[1] = dx = σ*(y-x)
    du[2] = dy = x*(ρ-z) - y
    du[3] = dz = x*y - β*z
end

u0 = [1.0, 0.0, 0.0]
tspan = (0.0, 100.0)
p = [10.0, 28.0, 8/3]

prob = ODEProblem(parametrized_lorenz!, u0, tspan, p)

sol = solve(prob)
plot(sol, idxs=(1,2,3))


# rich output for saving with nice variable names
f! = ODEFunction(parametrized_lorenz!,syms=[:x,:y,:z])
prob = ODEProblem(f!, u0, tspan, p)
sol = solve(prob)
df = DataFrame(sol)


## integrator interface
dt = 0.1
integrator = init(prob; u0, tspan, p, saveat=dt, save_on=false)

fieldnames(typeof(integrator))

n = 20
for i ∈ 1:n
    step!(integrator, dt, true)
    println(integrator.t)
    @assert check_error(integrator) == ReturnCode.Success
end

integrator.sol
integrator.t
integrator.u


