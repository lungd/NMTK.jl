using NMTK

@parameters t x
ps = collect(@parameters R D GL EL)
vars = collect(@variables v(..))
Dt = Differential(t)
Dxx = Differential(x)^2

eqs = [
    Dt(v(t,x)) ~ GL*(EL-v(t,x)) + (D / 4R)*(Dxx(v(t,x))) + 0.0
]

defaults = Dict(R => 40.0,
                D => 10.0,
                GL => 0.1,
                EL => -65.0)


bcs = [v(0.0,x) ~ -65.0,
       v(t,0.0) ~ 15.0,# for all t > 0
       #v(t,0.0) ~ 15.0,# for all t > 0
       ] #for all  0 < x < 1]

domains = [t ∈ IntervalDomain(0.0,100.0),
           x ∈ IntervalDomain(0.0,10.0)]

sys = PDESystem(eqs,bcs,domains,[t,x],[v(t,x)],ps,defaults)

dx = 0.1
order = 2
dx = range(0.0,10,length=100)
#@time prob = discretize(sys,discretization)

discretization = MOLFiniteDifference([x=>0.1],t)
prob = discretize(sys,discretization)


discretization = MOLFiniteDifference([x=>dx],t)
discretization_centered = MOLFiniteDifference([x=>dx],t;centered_order=order)
# Higher order centered difference
discretization_centered_order4 = MOLFiniteDifference([x=>dx],t;centered_order=4)

# Convert the PDE problem into an ODE problem
prob = discretize(sys,discretization)
prob_centered = discretize(sys,discretization_centered)
prob_centered_order4 = discretize(sys,discretization_centered_order4)

@time sol = solve(prob,Tsit5())
xs = collect(range(0.0,10,length=length(sol.u[1])))
ts = sol.t
vs = transpose(hcat(sol.u...))
heatmap(xs,ts,vs)

@time sol = solve(prob_centered,Tsit5())
xs = collect(range(0.0,10,length=length(sol.u[1])))
ts = sol.t
vs = transpose(hcat(sol.u...))
heatmap(xs,ts,vs)

@time sol = solve(prob_centered_order4,Tsit5())
xs = collect(range(0.0,10,length=length(sol.u[1])))
ts = sol.t
vs = transpose(hcat(sol.u...))
heatmap(xs,ts,vs)
