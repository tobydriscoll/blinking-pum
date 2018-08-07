% INPUT:
%      PUApprox: PUApprox approximation
%             v: given solution
%         evalF: residual function for leaf
%           Jac: jacobian function for leaf
%
% OUTPUT:
%          r_er: correction of solution
%          J_v_pls_er: jacobian of J_hat(v+er)
%
% NOTE sol is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then sol = [u1;u2;v1;v2].
function [ r_er,J_v_pls_er ] = CoarseCorrect( PUApprox, v,evalf,jacf)


num_sols = length(v)/length(PUApprox);

v_hat = [];
for i=1:num_sols
    v_hat = [v_hat;PUApprox.Fine2Coarse(v((1:length(PUApprox)) + (i-1)*length(PUApprox)))];
end

r = ParResidual(v,PUApprox,evalf);

r_hat = [];
for i=1:num_sols
    r_hat = [r_hat;PUApprox.Fine2Coarse(r((1:length(PUApprox)) + (i-1)*length(PUApprox)))];
end

PUApprox.Coarsen();

r_hat = r_hat - ParResidualFun(v_hat,PUApprox,evalf);

params = [20,-1,.5,0];
tol = [1e-10 1e-11];

[er,~,~,~,~] = nsoldAS(rand(size(v_hat)),@(er)Residual(er,v_hat,r_hat,PUApprox,evalf),@(er)CoarseASJac(PUApprox,jacf,er,v_hat),tol,params);

J_v_pls_er = CoarseASJac(PUApprox,jacf,er,v_hat);

r_er = [];

for i=1:num_sols
    r_er = [r_er;PUApprox.Coarse2Fine(er((1:length(PUApprox)) + (i-1)*length(PUApprox)))];
end

PUApprox.Refine();

end

function F = Residual(er,v,r,PUApprox,evalf)
    F = ParResidualFun(v+er,PUApprox,evalf);
    F = ParResidualInterface(F,er,PUApprox);
    F = F + r;
end

