
using JuMP
using PiecewiseLinearOpt

# Define data
include("data.jl")
@assert length(h_ρi) == length(ρi) == length(Hi)
@assert issorted(h_ρi)
@assert issorted(ρi, rev=true)
@assert issorted(Hi)

# Define global model in log-space
function global_model(solver, n_pts::Int)
    m = Model(solver=solver)

    # log variables
    @variables m begin
        # Subsystem specific parameters
        D_p # payload optical aperture diameter
        D_T # transmitter antenna diameter
        A # solar panel surface area
        ρ_A
        η_A

        # Mass related parameters
        m_t >= m_c # total mass
        m_b # battery mass
        m_A # solar panel mass
        m_p # payload mass
        m_T # transmitter mass
        m_S # structural mass
        m_P # propulsion mass

        # Power related parameters
        P_t # total power
        P_T # transmitter power
        E_b # battery energy

        # Orbit related parameters
        h_min <= h <= h_max # altitude
        a # semi-major axis
        T # orbit period
        r # worst case communications distance
        d_min <= d <= d_max # daylight fraction of orbit
        #e_min <= e <= e_max # eclipse fraction of orbit
        g_min <= g <= g_max # ground station viewing fraction of orbit

        h_ρi[1] <= h_ρ <= h_ρi[end]

        # Lifetime (years)
        exp_Ln_min <= exp_Ln <= exp_Ln_max # without propulsion (orbit decay)
        exp_Lp_min <= exp_Lp <= exp_Lp_max # with propulsion
        ρ # atmospheric density
        H # atmosphere scale height

        # Auxiliary variables
        exp_pTt
        exp_Plt
        exp_Ra
        exp_ha
        exp_hr
        exp_Rhr
        exp_mbt
        exp_mpt
        exp_mAt
        exp_mTt
        exp_mPt
        exp_mSt
        exp_mct
        exp_d

        # Component related variables
        x_T[1:n_T], Bin
        x_b[1:n_b], Bin
        x_p[1:n_p], Bin
        x_A[1:n_A], Bin
    end

    # Catalog selections
    @constraints m begin
        sum(x_T) == 1 # transmitter catalog
        D_T == dot(x_T, D_Ti)
        m_T == dot(x_T, m_Ti)
        P_T == dot(x_T, P_Ti)
        #G_T == dot(x_T, G_Ti)

        sum(x_b) == 1 # battery catalog
        E_b == dot(x_b, E_bi)
        m_b == dot(x_b, m_bi)

        sum(x_p) == 1 # payload catalog
        D_p == dot(x_p, D_pi)
        m_p == dot(x_p, m_pi)

        sum(x_A) == 1 # solar panel catalog
        ρ_A == dot(x_A, ρ_Ai)
        η_A == dot(x_A, η_Ai)
    end

    # Piecewiselinear function values from data table
    ρ = piecewiselinear(m, h_ρ, h_ρi, ρi)
    H = piecewiselinear(m, h_ρ, h_ρi, Hi)

    # Piecewiselinear function approximations
    brk_log = (v_min, v_max, num_pts) -> log.(linspace(exp(v_min), exp(v_max), num_pts))
    # Orbit
    pwgraph_fhR = piecewiselinear(m, h, brk_log(h_min, h_max, n_pts), h_val -> log(acos(1/(exp(h_val - R) + 1)))) # concave
    pwgraph_exp_d = piecewiselinear(m, d, brk_log(d_min, d_max, n_pts), exp) # convex
    #pwgraph_exp_e = piecewiselinear(m, e, brk_log(e_min, e_max, n_pts), exp) # convex
    pwgraph_exp_g = piecewiselinear(m, g, brk_log(g_min, g_max, n_pts), exp) # convex
    # Lifetime
    pwgraph_Ln = piecewiselinear(m, exp_Ln, linspace(exp_Ln_min, exp_Ln_max, n_pts), log) # concave
    pwgraph_Lp = piecewiselinear(m, exp_Lp, linspace(exp_Lp_min, exp_Lp_max, n_pts), log) # concave

    # Convex extended formulation constraints
    @NLconstraints m begin
        exp(P_T - P_t) <= exp_pTt
        exp(P_l - P_t) <= exp_Plt
        exp(R - a) <= exp_Ra
        exp(h - a) <= exp_ha
        exp(2*h - 2*r) <= exp_hr
        exp(log(2) + R + h - 2*r) <= exp_Rhr
        # Mass
        exp(m_b - m_t) <= exp_mbt
        exp(m_p - m_t) <= exp_mpt
        exp(m_A - m_t) <= exp_mAt
        exp(m_T - m_t) <= exp_mTt
        exp(m_P - m_t) <= exp_mPt
        exp(m_S - m_t) <= exp_mSt
        exp(m_c - m_t) <= exp_mct
        # Orbit
        exp(d) <= exp_d
        log(π) + g <= log(acos(1/(exp(h - R) + 1))) # instead of linearized
    end

    # Linear constraints
    @constraints m begin
        # Power and communications
        P_t - d == A + η_A + Q
        E_b >= P_t - d + T
        EN + L + k + T_s + R + log(2π) + N + B + log(4) + 2*r == P_T + G_r + X_r + g + T + η + 2*D_T
        # Payload performance
        X_r == h + λ_v - D_p
        # Orbit
        # g == α_1 + γ_1*(h - R) # linearized
        T == log(2π) + (3*a - μ)/2
        # Mass budgets
        #m_b == ρ_b + E_b
        #m_p == ρ_p + 3/2*D_p
        m_A == ρ_A + A
        #m_T == ρ_T + 3/2*D_T
        #m_P == ρ_P - h
        m_S == η_S + m_t
        # Lifetime
        pwgraph_Ln + log(3600*24*365) == H + m_t + T - log(2π) - C_D - A - 2*h - ρ
        pwgraph_Lp + log(3600*24*365) == m_P + I_sp + G + a - log(0.5) - C_D - A - ρ - μ
        exp_Ln + exp_Lp >= exp_Lt_min
        # From convex constraints
        exp_pTt + exp_Plt <= 1
        exp_Ra + exp_ha <= 1
        exp_hr + exp_Rhr <= 1
        exp_mbt + exp_mpt + exp_mAt + exp_mTt + exp_mPt + exp_mSt + exp_mct <= 1
        # From nonconvex constraints
        exp_d <= pwgraph_exp_d
        #pwgraph_exp_e + exp_d >= 1
        exp_d <= pwgraph_exp_g + 1/2
        log(π) + g >= pwgraph_fhR
    end

    # Minimize total mass
    @objective(m, Min, m_t)

    # Solve
    status = solve(m)
    println()
    println("  status: ", status)
    println("  total mass: ", exp(getobjectivevalue(m)))
    println()
    println("  d = ", exp(getvalue(d)))
    #println("  e = ", exp(getvalue(e)))
    println("  g = ", exp(getvalue(g)))
    println()
    println("  Ln = ", getvalue(exp_Ln))
    println("  Lp = ", getvalue(exp_Lp))
    println()
    println("  transmitter: ", find(i -> (i > 0.5), getvalue(x_T)))
    println("  battery: ", find(i -> (i > 0.5), getvalue(x_b)))
    println("  payload: ", find(i -> (i > 0.5), getvalue(x_p)))
    println("  solar panel: ", find(i -> (i > 0.5), getvalue(x_A)))
    println()
    ;
end

# Define solvers
using CPLEX
using MINLPOA
global_solver = MINLPOASolver(log_level=1, mip_solver=CplexSolver(CPX_PARAM_SCRIND=1, CPX_PARAM_EPINT=1e-9, CPX_PARAM_EPRHS=1e-9, CPX_PARAM_EPGAP=1e-7))

# Run
global_model(global_solver, 1000)



# using Ipopt
# using Pajarito
# global_solver = PajaritoSolver(
#     mip_solver=CplexSolver(CPX_PARAM_SCRIND=1, CPX_PARAM_EPINT=1e-9, CPX_PARAM_EPRHS=1e-9, CPX_PARAM_EPGAP=1e-7),
#     cont_solver=IpoptSolver(print_level=0),
#     mip_solver_drives=true,
#     log_level=1,
#     rel_gap=1e-7)
