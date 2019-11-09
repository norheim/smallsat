using JuMP
using Ipopt
solver = IpoptSolver(print_level=5)
m = Model(solver=solver)
n = 3
@variables m begin
    x[1:n, 1:n], Bin
    y[1:n], Bin
    z[1:n, 1:n], Bin
end
@constraints m begin
    for i=1:n
        sum(x[i]) == 1
    end
    for i=1:n
        sum(z[i]) == 1
    end
    y[1]-y[2]<=n*x[1, 1]
    y[1]-y[2]<=n*x[1, 1]
end

status = solve(m)
