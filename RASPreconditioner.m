function [ output] = RASPreconditioner(A,b,PUApprox,domain,force,border,sol )
Residual = b - A(sol);
output = sol + RASStep(PUApprox,domain,force,border,Residual);
end