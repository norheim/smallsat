
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
        -2 <= P_T <= 2 # transmitter power
        exp_pTt >= 0
        x_T[1:3], Bin
    end

    # Catalog selections
    @constraints m begin
        sum(x_T) == 1 # transmitter catalog
        P_T == dot(x_T, [-1,0,1])
    end

    # Convex extended formulation constraints
    @NLconstraints m begin
        -exp(P_T) >= -exp_pTt
        exp(P_T) <= exp_pTt
    end

    @objective(m, Min, exp_pTt)

    # Solve
    status = solve(m)
    println()
    println("  status: ", status)
    println("  total: ", exp(getobjectivevalue(m)))
    ;
end

# Define solvers
using Ipopt
using CPLEX
using Pajarito
using MINLPOA

# global_solver = PajaritoSolver(
#     mip_solver=CplexSolver(CPX_PARAM_SCRIND=1, CPX_PARAM_EPINT=1e-9, CPX_PARAM_EPRHS=1e-9, CPX_PARAM_EPGAP=1e-7),
#     cont_solver=IpoptSolver(print_level=0),
#     mip_solver_drives=true,
#     log_level=1,
#     rel_gap=1e-7)

global_solver = MINLPOASolver(log_level=1, mip_solver=CplexSolver(CPX_PARAM_SCRIND=1, CPX_PARAM_EPINT=1e-9, CPX_PARAM_EPRHS=1e-9, CPX_PARAM_EPGAP=1e-7))

# Run
global_model(global_solver, 60)
