% CoarseCorrect
% This method solves for the nonlinear coarse correction used in the 2
% level SNK method.
% 
% INPUT:
%      PUApprox: PUApprox approximation
%             v: given solution
%
%             f: nonlinear residual function f(x,p) for solution x and local
%             approximation p. f(x,p) evaluates the residual on the domain
%             of p.
%
%            Jac: Jacobian function Jac(x,p) for solution x and local
%             approximation p. Jac(x,p) evaluates the residual on the domain
%             of p.
%
% OUTPUT:
%          r_er: correction of solution
%          J_v_pls_er: jacobian of J_hat(v+er)
%
% NOTE sol is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then sol = [u1;u2;v1;v2].
function [ r_er,J_v_pls_er ] = CoarseCorrect( PUApprox, v,f,Jac,j,tol_c)


num_sols = length(v)/length(PUApprox);

v_hat = [];
for i=1:num_sols
    v_hat = [v_hat;PUApprox.Fine2Coarse(v((1:length(PUApprox)) + (i-1)*length(PUApprox)),j)];
end

r = ParResidual(v,PUApprox,f);

r_hat = [];
for i=1:num_sols
    r_hat = [r_hat;PUApprox.Fine2Coarse(r((1:length(PUApprox)) + (i-1)*length(PUApprox)),j)];
end

PUApprox.Coarsen();

r_hat = r_hat - ParResidual(v_hat,PUApprox,f);

RES = @(u) ParResidual(u,PUApprox,f)+r_hat;
JAC = @(u)CoarseASJac(PUApprox,Jac,u);

%options = optimoptions(@fsolve,'SpecifyObjectiveGradient',true,'MaxIterations',50,'FunctionTolerance',tol_c,'Display','off');
options = optimoptions(@fsolve,'SpecifyObjectiveGradient',true,'MaxIterations',50,'Display','off');
er = fsolve(@(u)sol_and_jac(@(u)RES(u),@(u)JAC(u),u),v_hat,options);

J_v_pls_er = JAC(er);

er = er(:,end) - v_hat;



r_er = [];

for i=1:num_sols
    r_er = [r_er;PUApprox.Coarse2Fine(er((1:length(PUApprox)) + (i-1)*length(PUApprox)))];
end

PUApprox.Refine();

end

function F = Residual(er,v,r,PUApprox,evalf)
    F = ParResidual(v+er,PUApprox,evalf);
    F = F + r;
end

