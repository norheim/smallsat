# include("C:\\Users\\johan\\PycharmProjects\\smallsat\\gsstation.jl")
using JuMP, PiecewiseLinearOpt, Ipopt, Gurobi, Pavito
mip_solver = Gurobi.GurobiSolver(OutputFlag=1, IntFeasTol=1e-9, 
    FeasibilityTol=1e-9, MIPGap=1e-7)

global_solver = PavitoSolver(
    mip_solver=mip_solver,
    cont_solver=IpoptSolver(print_level=0),
    mip_solver_drives=true,
    log_level=1,
    rel_gap=1e-7,
    )

R = 6378e3 # m
μ = 3.986005e14
h = 400e3
a = h+R
K = log(86400^2*μ/(4*π^2))
T2 = 24*3600
T1 = 2π*sqrt(a^3/μ)
loga = log(a)

m = Model(solver=global_solver)

orbits = 1
N = 20*orbits
n = 1
maxn = 1
srand(100)
lat_g, lon_g =  (0.5-rand(n))*π/2, rand(n)*2π
t_exp = log.(linspace(0.1,540*N, N)) # one step per 60 sec
x = linspace(0,2π,N)*orbits # one full satellite orbit

inc = -0.97
sun_cd = sqrt.(sin.(x).^2+sin(inc)^2*cos.(x).^2)
sun_cd2 = Int.(cos.(x) .>= 0)

knum = T2/T1
h_min = 300e3
h_max = 500e3
a_min = h_min+R
a_max = h_max+R
k_max = 24*3600/(2π*sqrt(a_min^3/μ))
k_min = 24*3600/(2π*sqrt(a_max^3/μ))

brk = (x_min, x_max, num_pts) -> 
    log.(linspace(exp(x_min), exp(x_max), num_pts))
brk_log = (x_min, x_max, num_pts) -> 
    log.(linspace(x_min, x_max, num_pts))

#@variable(m, -π/2<=inc<=π/2)
#@variable(m, k)
@variable(m, k_exp)
@variable(m, log(h_min)<=h_exp<=log(h_max))
@variable(m, a_min<=a<=a_max)
@variable(m, log(a_min)<=a_exp<=log(a_max))
@variable(m, x_v[1:N, 1:n], Bin)
@variable(m, x_s[1:N], Bin)
@variable(m, x_i[1:n], Bin)

#sinc = piecewiselinear(m, inc, -π/2:0.2:π/2, (inc) -> sin(inc))
#cinc = piecewiselinear(m, inc, -π/2:0.2:π/2, (inc) -> cos(inc))


npts = 10
a_exp = piecewiselinear(m, h_exp, brk_log(h_min, h_max, npts), 
       (h)->log(exp(h) + R))
a = piecewiselinear(m, a_exp, brk_log(a_min, a_max, npts), 
            (a)->exp(a))
    
d = [piecewiselinear(m, k_exp, a_exp, brk_log(k_min,k_max,npts), 
        brk_log(a_min, a_max, npts),
        (kexp, a) -> exp(a)/R*(-sin(inc)*sin(lat_g[j])*cos(x[i])
        +sin(x[i])*sin(lon_g[j]+x[i]/exp(kexp))*cos(lat_g[j])
        +cos(inc)*cos(lat_g[j])*cos(x[i])*cos(lon_g[j]+x[i]/exp(kexp))))-1
    for i=1:N for j=1:n]

@objective(m, Max, sum(x_v))
@constraint(m, [i=1:N,j=1:n],  d[i,j] <= 4*x_v[i,j])
@constraint(m, [i=1:N,j=1:n],  d[i,j] >= -4*(1-x_v[i,j]))
@constraint(m, [i=1:N], a/R*(sun_cd[i]+sun_cd2[i]) - 1 <= 4*x_s[i] )
@constraint(m, [i=1:N], a/R*(sun_cd[i]+sun_cd2[i]) - 1>= -4*(1-x_s[i]))

#@NLconstraint(m, exp(a) <= a_exp)
#@constraint(m, k_exp == log(knum))
@constraint(m, log(2π)+k_exp+1/2*(3*a_exp-log(μ)) == log(T2))

@constraint(m, [i=1:N, j=1:n], x_v[i,j] <= x_i[j]) # which ground stations we select
@constraint(m, sum(x_i) <= maxn)

# status = solve(m)