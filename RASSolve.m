function [ output ] = RASSolve(PUApprox,domain,force,border,sol)
is_converged = false;
tol = 1e-5;
while ~is_converged
    new_sol = RASPreconditioner(PUApprox,domain,force,border,sol);
    is_converged = norm(new_sol-sol,inf)/norm(sol,inf)<tol;
    sol = new_sol;
end
output = sol;
end

