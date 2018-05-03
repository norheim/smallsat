
using JuMP

# Define data
include("data.jl")

# Define local model in log-space
function local_model(solver)
    m = Model(solver=solver)

    # log variables
    @variables m begin
        # Subsystem specific parameters
        D_p # payload optical aperture diameter
        D_T # transmitter antenna diameter
        A # solar panel surface area

        # Mass related parameters
        m_t # total mass
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
        h_min <= h <= h_max, (start=log(1.0386e6)) # altitude
        a # semi-major axis
        T # orbit period
        r # worst case communications distance
        d_min <= d <= d_max, (start=log(0.6723)) # daylight fraction of orbit
        e_min <= e <= e_max, (start=log(0.3277)) # eclipse fraction of orbit
        g_min <= g <= g_max, (start=log(0.1723)) # ground station viewing fraction of orbit
    end

    # Nonconvex constraints
    @NLconstraints m begin
        # Orbit
        exp(e) >= 1 - exp(d)
        exp(g) >= exp(d) - 1/2
        log(π) + g == log(acos(1/(exp(h - R) + 1))) # instead of linearized
    end

    # Convex constraints
    @NLconstraints m begin
        # Power and communications
        exp(P_T - P_t) + exp(P_l - P_t) <= 1
        # Orbit
        exp(R - a) + exp(h - a) <= 1
        exp(2*h - 2*r) + exp(log(2) + R + h - 2*r) <= 1
        # Mass budgets
        exp(m_b - m_t) + exp(m_p - m_t) + exp(m_A - m_t) + exp(m_T - m_t) + exp(m_P - m_t) + exp(m_S - m_t) + exp(m_c - m_t) <= 1
    end

    # Linear constraints
    @constraints m begin
        # Power and communications
        P_t - d == A + η_A + Q
        E_b >= P_t - d + e + T
        EN + L + k + T_s + R + log(2π) + N + B + log(4) + 2*r == P_T + G_r + X_r + g + T + η + 2*D_T
        # Payload performance
        X_r == h + λ_v - D_p
        # Orbit
        # g == α_1 + γ_1*(h - R) # linearized
        T == log(2π) + (3*a - μ)/2
        # Mass budgets
        m_b == ρ_b + E_b
        m_p == ρ_p + 3/2*D_p
        m_A == ρ_A + A
        m_T == ρ_T + 3/2*D_T
        m_P == ρ_P - h
        m_S == η_S + m_t
    end

    # Minimize total mass
    @objective(m, Min, m_t)

    # Solve
    status = solve(m)
    println("\nlocal:")
    println("  status: ", status)
    println("  total mass: ", exp(getobjectivevalue(m)))
    println("  d = ", exp(getvalue(d)))
    println("  e = ", exp(getvalue(e)))
    println("  g = ", exp(getvalue(g)))
    println()

    # @show exp(getvalue(D_p))
    # @show exp(getvalue(D_T))
    # @show exp(getvalue(A))
    # @show exp(getvalue(m_t))
    # @show exp(getvalue(m_b))
    # @show exp(getvalue(m_A))
    # @show exp(getvalue(m_p))
    # @show exp(getvalue(m_T))
    # @show exp(getvalue(m_S))
    # @show exp(getvalue(m_P))
    # @show exp(getvalue(P_t))
    # @show exp(getvalue(P_T))
    # @show exp(getvalue(E_b))
    # @show exp(getvalue(h))
    # @show exp(getvalue(a))
    # @show exp(getvalue(T))
    # @show exp(getvalue(r))
    # @show exp(getvalue(d))
    # @show exp(getvalue(e))
    # @show exp(getvalue(g))
    ;
end

# Define solvers
using Ipopt
local_solver = IpoptSolver(print_level=5)

# Run
local_model(local_solver)
