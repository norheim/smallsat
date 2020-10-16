# include("C:\\Users\\johan\\PycharmProjects\\smallsat\\ordesign.jl")
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

# Define data
include("data.jl")

X_R = 20
B = 32
n_pts = 10
exp_Lt_min = 3

X_R_original = X_R
X_R = log(X_R)
B = log(B)
m = Model(solver=global_solver)

# log variables
@variables m begin
    # Subsystem specific parameters
    D_p # payload optical aperture diameter
    X_r # payload resolution
    G_T # transmitter antenna diameter
    A # solar panel surface area

    # catalog
    ρ_A
    η_A

    # Mass related parameters
    #X_r
    m_t >= m_min # total mass
    m_b # battery mass
    m_A # solar panel mass
    m_p # payload mass
    m_T # transmitter mass
    m_S # structural mass
    m_P # propulsion mass

    # Comms related parameters
    dRate
    data
    fr_comms
    fr_charge
    expfrcm
    expfrch
    

    # Power related parameters
    P_p # payload power
    P_t # total power
    P_T # transmitter power
    E_b # battery energy

    # Orbit related parameters
    h_min <= h <= h_max # altitude
    a # semi-major axis
    T # orbit period
    r # worst case communications distance
    d_min <= d <= d_max # daylight fraction of orbit
    e_min <= e <= e_max # eclipse fraction of orbit
    g_min <= g <= g_max # ground station viewing fraction of orbit

    # Lifetime (years)
    exp_Ln_min <= exp_Ln <= exp_Ln_max # without propulsion (orbit decay)
    exp_Lp_min <= exp_Lp <= exp_Lp_max # with propulsion
    ρ # atmospheric density
    H # atmosphere scale height
    h_ρi[1] <= h_ρ <= h_ρi[end]

    # Disjunctive
    T_g
    log(1e-2) <= m_P2 <= log(10) # mass for reaction wheel
    log(1e-2) <= m_M <= log(10) # mass for magnetorques

    # Auxiliary variables
    exp_PTt
    exp_Plt
    #exp_Ra
    #exp_ha
    #exp_hr
    #exp_Rhr
    exp_mbt
    exp_mpt
    exp_mAt
    exp_mTt
    exp_mPt
    exp_mSt
    exp_mct
    exp_d
    a_exp

    # Disjunctive
    exp_mMPt >= 0.0

    # Component related variables
    #x_P2, Bin
    #x_M, Bin
    x_T[1:n_T], Bin
    x_b[1:n_b], Bin
    x_p[1:n_p], Bin
    x_A[1:n_A], Bin

    # Operations
    k_exp
    x_v[1:Nt, 1:n_gs], Bin
    x_vc[1:Nt, 1:n_gs], Bin
    x_s[1:Nt], Bin
    x_sc[1:Nt], Bin
    x_i[1:Nt], Bin
    x_pc[1:Nt], Bin
end

# Catalog selections
@constraints m begin
    #x_P2 + x_M == 1

    sum(x_T) == 1 # transmitter catalog
    G_T == dot(x_T, G_Ti)
    m_T == dot(x_T, m_Ti)
    P_T == dot(x_T, P_Ti)
    #G_T == dot(x_T, G_Ti)
    #m_T == ρ_T + 3/2*D_T

    sum(x_b) == 1 # battery catalog
    E_b == dot(x_b, E_bi)
    m_b == dot(x_b, m_bi)
    #m_b == ρ_b + E_b

    sum(x_p) == 1 # payload catalog
    D_p == dot(x_p, D_pi)
    m_p == dot(x_p, m_pi)
    P_p == dot(x_p, P_pi)
    #m_p == ρ_p + 3/2*D_p

    sum(x_A) == 1 # solar panel catalog
    ρ_A == dot(x_A, ρ_Ai)
    η_A == dot(x_A, η_Ai)
end

# Piecewiselinear function values from data table
ρ = piecewiselinear(m, h_ρ, h_ρi, ρi)
H = piecewiselinear(m, h_ρ, h_ρi, Hi)

# Nonconvex extended formulation piecewiselinear approx
brk_log = (v_min, v_max, num_pts) -> log.(linspace(exp(v_min), exp(v_max), num_pts))
brk = (x_min, x_max, num_pts) -> log.(linspace(x_min, x_max, num_pts))
pwgraph_exp_d = piecewiselinear(m, d, brk_log(d_min, d_max, n_pts), exp) # convex
pwgraph_exp_e = piecewiselinear(m, e, brk_log(e_min, e_max, n_pts), exp) # convex
pwgraph_exp_g = piecewiselinear(m, g, brk_log(g_min, g_max, n_pts), exp) # convex
fhR = h_val -> log(acos(1/(exp(h_val - R) + 1))) # instead of linearized
pwgraph_fhR = piecewiselinear(m, h, brk_log(h_min, h_max, n_pts), fhR) # concave
ahr = h_val -> log(exp(h_val) + exp(R))
pwgraph_a = piecewiselinear(m, h, brk_log(h_min, h_max, n_pts), ahr) #convex
a_exp = piecewiselinear(m, h, brk_log(h_min, h_max, n_pts), 
    (h_val)-> exp(h_val) + exp(R))

# Lifetime
pwgraph_Ln = piecewiselinear(m, exp_Ln, linspace(exp_Ln_min, exp_Ln_max, n_pts), log) # concave
pwgraph_Lp = piecewiselinear(m, exp_Lp, linspace(exp_Lp_min, exp_Lp_max, n_pts), log) # concave

# Operations
θ = linspace(0, 2π, Nt)*orbits # one full satellite orbit

inc = -0.97
sun_cd = sqrt.(sin.(θ).^2+sin(inc)^2*cos.(θ).^2)
sun_cd2 = Int.(cos.(θ) .>= 0)

gs_cd = [piecewiselinear(m, k_exp, pwgraph_a, brk(k_min,k_max,n_pts), 
        brk_log(a_min, a_max, n_pts),
        (kexp, a) -> exp(a)/exp(R)*(-sin(inc)*sin(lat_g[j])*cos(θ[i])
        +sin(θ[i])*sin(lon_g[j]+θ[i]/exp(kexp))*cos(lat_g[j])
        +cos(inc)*cos(lat_g[j])*cos(θ[i])*cos(lon_g[j]+θ[i]/exp(kexp))))-1
    for i=1:Nt for j=1:n_gs]


@constraint(m, k_exp + T == log(T2))
@constraint(m, [i=1:Nt,j=1:n_gs],  gs_cd[i,j] <= 4*x_v[i,j])
@constraint(m, [i=1:Nt,j=1:n_gs],  gs_cd[i,j] >= -4*(1-x_v[i,j]))
@constraint(m, [i=1:Nt], a_exp/exp(R)*(sun_cd[i]+sun_cd2[i]) - 1 <= 4*x_s[i] )
@constraint(m, [i=1:Nt], a_exp/exp(R)*(sun_cd[i]+sun_cd2[i]) - 1 >= -4*(1-x_s[i]))
@constraint(m, [i=1:Nt, j=1:n_gs], x_v[i,j] <= x_i[j]) # which ground stations we select
@constraint(m, sum(x_i) <= max_n_gs)
@constraint(m, x_sc .<= x_s)
@constraint(m, x_vc .<= x_v)
@constraint(m, [i=1:Nt], x_sc[i]+sum(x_vc[i,:]) <= 1) # cannot do sun and comms
@constraint(m, [i=1:Nt], sum(x_vc[i,:])+x_pc[i] <= 1) # cannot do payload and comms

expfrcm = piecewiselinear(m, fr_comms, linspace(0.9*1/Nt, 0.9, n_pts), 
    (fr_comms)-> log(fr_comms))
expfrch = piecewiselinear(m, fr_charge, linspace(0.9*1/Nt, 0.9, n_pts), 
    (fr_charge)-> log(fr_charge))

r = piecewiselinear(m, h, brk_log(h_min, h_max, n_pts), 
        (h) -> log(sqrt(exp(h)^2+2*exp(R)*exp(h))))

@constraints m begin
    Nt*fr_comms == sum(x_vc)
    dRate == data - expfrcm - T

    Nt*fr_charge == sum(x_sc)
    #fr_charge >= 0.5
    #A >= log(0.1)
    P_t - expfrch == A + η_A + Q
    #fr_payload == sum(x_p)/20
    #data_collected = data + expfrp + T
    #expPc + expPp <= 1
end

# Convex extended formulation constraints
@NLconstraints m begin
    #exp(a) <= a_exp
    #exp(2*h - 2*r) <= exp_hr
    #exp(log(2) + R + h - 2*r) <= exp_Rhr

    exp(d) <= exp_d
    log(π) + g <= log(acos(1/(exp(h - R) + 1)))

    #exp(m_P2 - m_t) <= exp_mMPt
    exp(m_M - m_t) <= exp_mMPt
    #(x_P2 + 1e-3)*exp((m_P2 - m_t)/(x_P2 + 1e-3)) <= exp_mMPt
    #(x_M + 1e-3)*exp((m_M - m_t)/(x_M + 1e-3)) <= exp_mMPt
    
    #exp(fr_comms) <= expfrcm
    #exp(fr_charge) <= expfrch
    #exp(fr_payload) <= expfrp
    #exp(P_T - P_t + expfrcm - expfrch) <= expPc
    #exp(P_T - P_p + expfrcm - expfrp) <= expPp
    
    exp(P_T - P_t) <= exp_PTt
    exp(P_l - P_t) <= exp_Plt

    exp(m_b - m_t) <= exp_mbt
    exp(m_p - m_t) <= exp_mpt
    exp(m_A - m_t) <= exp_mAt
    exp(m_T - m_t) <= exp_mTt
    exp(m_P - m_t) <= exp_mPt
    exp(m_S - m_t) <= exp_mSt
    exp(m_c - m_t) <= exp_mct
end

# Linear constraints
@constraints m begin
    # Power and communications
    #P_t - d == A + η_A + Q
    E_b >= P_t - d + T
    P_T <= P_t
    data == log(2π) + R + N + B - X_r
    #dRate == data - g - T
    EN + L + k + T_s + dRate + log(4) + 2*r <= P_T + G_r + G_T + 2*(λ_c-log(π))
    # Payload performance
    X_R == X_r
    X_r == log(1.22) + h + λ_v - D_p
    # Orbit
    # g == α_1 + γ_1*(h - R) # linearized
    a == pwgraph_a
    T - log(2π) == (3*a - μ)/2

    # Mass budgets
    m_A == ρ_A + A
    #m_P == ρ_P - h
    m_S == η_S + m_t

    # lifetime
    pwgraph_Ln + log(3600*24*365) == T + H + m_t - log(2π) - C_D - A - 2*a - ρ
    pwgraph_Lp + log(3600*24*365) == m_P + I_sp + G + a - log(0.5) - C_D - A - ρ - μ
    exp_Ln + exp_Lp >= exp_Lt_min
    h_ρ == h

    # From convex constraints
    exp_PTt + exp_Plt <= 1
    #exp_Ra + exp_ha <= 1 # not tight!
    #exp_hr + exp_Rhr <= 1
    exp_mbt + exp_mpt + exp_mPt + exp_mAt + exp_mTt + exp_mSt + exp_mct + exp_mMPt<= 1

    # From nonconvex constraints
    exp_d <= pwgraph_exp_d
    pwgraph_exp_e + exp_d == 1
    exp_d <= pwgraph_exp_g + 1/2
    log(π) + g >= pwgraph_fhR

    # Momentum budgets
    T_g == log(3) + μ + c_W - log(2) - 3*a
    #m_P2 == ρ_P2 + T + T_g + log(1/4*sqrt(2)/2*3*365) + log(exp_Lt_min) - I_sp - G
    m_M == ρ_M + T_g
end

# Minimize total mass
@objective(m, Min, m_t)

# Solve
status = solve(m)
D = [
    ("status", status),
    ("m_t" , exp(getvalue(m_t))),
    ("x_p" , find(i -> (i > 0.5), getvalue(x_p))),
    ("x_R" , exp(X_R)),
    ("x_r" , exp(getvalue(X_r))),
    ("x_B" , find(i -> (i > 0.5), getvalue(x_b))),
    ("m_B" , exp(getvalue(m_b))),
    ("m_S" , exp(getvalue(m_S))),
    ("x_T" , find(i -> (i > 0.5), getvalue(x_T))),
    ("m_T" , exp(getvalue(m_T))),
    ("G_T" , exp(getvalue(G_T))),
    ("data" , 2*π*exp(R)/exp(getvalue(X_r))*exp(B)*exp(N)),
    ("rate" , 2*π*exp(R)/exp(getvalue(X_r))*exp(B)*exp(N)/exp(getvalue(T))),
    ("x_A" , find(i -> (i > 0.5), getvalue(x_A))),
    ("A" , exp(getvalue(A))),
    ("m_A" , exp(getvalue(m_A))),
    ("Ln" , getvalue(exp_Ln)),
    ("Lp" , getvalue(exp_Lp)),
    ("Lt" , getvalue(exp_Ln)+getvalue(exp_Lp)),
    ("h" , exp(getvalue(h))/1000),
    ("a" , exp(getvalue(a))/1000),
    ("r" , exp(getvalue(r))/1000),
    ("d" , exp(getvalue(d))),
    ("e" , exp(getvalue(e))),
    ("g" , exp(getvalue(g))),
    ("T" , exp(getvalue(T))/60),
    ("ρ" , exp(getvalue(ρ))),
    ("H" , exp(getvalue(H))),
    ("m_P" , exp(getvalue(m_P))),
    #"x_disj" , getvalue(x_P2) > 0.5 ? "P" : "M",
    #"x_P2" , getvalue(x_P2),
    #"x_M" , getvalue(x_M),
    ("T_g" , exp(getvalue(T_g))),
    #"m_P2" , exp(getvalue(m_P2)),
    ("m_M" , exp(getvalue(m_M))),
    ("P_t" , exp(getvalue(P_t))),
    ("P_T" , exp(getvalue(P_T))),
    ("exp_mbt" , getvalue(exp_mbt)),
    ("exp_mpt" , getvalue(exp_mpt)),
    ("exp_mPt" , getvalue(exp_mPt)),
    ("exp_mAt" , getvalue(exp_mAt)),
    ("exp_mTt" , getvalue(exp_mTt)),
    ("exp_mSt" , getvalue(exp_mSt)),
    ("exp_mct" , getvalue(exp_mct)),
    ("exp_mMPt" , getvalue(exp_mMPt)),
    #("exp_hr" , getvalue(exp_hr)),
    #("exp_Rhr" , getvalue(exp_Rhr)),
    ("mass constraint" , (getvalue(exp_mbt) + getvalue(exp_mpt) + getvalue(exp_mAt)
        + getvalue(exp_mTt) + getvalue(exp_mSt) + getvalue(exp_mct) + getvalue(exp_mMPt)+getvalue(exp_mPt))),
]