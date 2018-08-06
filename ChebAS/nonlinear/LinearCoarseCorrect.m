% INPUT:
%      PUApprox: PUApprox approximation        
%            w: solution to be applied to Jacobian
%               patches of f(v+e)
%      Jac_hat: f_hat'(v+e)-T_hat
%        T_hat: T_hat
%          FJv: cell array of f'(v_j)
%       FJ_hat: sparse matrix of f_hat'(v_hat)
%
% NOTE x is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then x = [u1;u2;v1;v2].
function [ y ] = LinearCoarseCorrect( PUApprox, w,Jac_hat,T_hat,FJv,FJv_hat)


num_sols = length(w)/length(PUApprox);

w_hat = [];
for i=1:num_sols
    w_hat = [w_hat;PUApprox.Fine2Coarse(w((1:length(PUApprox)) + (i-1)*length(PUApprox)))];
end

%need to linearize this
r = ParLinearResidual(w,PUApprox,FJv);

r_hat = [];

for i=1:num_sols
    r_hat = [r_hat;PUApprox.Fine2Coarse(r((1:length(PUApprox)) + (i-1)*length(PUApprox)))];
end

PUApprox.Coarsen();

b_hat = (Jac_hat+T_hat-FJv_hat)*w_hat+r_hat;

y_hat = -(Jac_hat)\b_hat;

y = [];

for i=1:num_sols
    y = [y;PUApprox.Coarse2Fine(y_hat((1:length(PUApprox)) + (i-1)*length(PUApprox)))];
end

PUApprox.Refine();

end
