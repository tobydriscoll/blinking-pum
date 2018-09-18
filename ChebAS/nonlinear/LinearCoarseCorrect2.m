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
function [ r_er ] = LinearCoarseCorrect2( PUApprox,w,g_mat,g_p_mat,J_er,T_hat,j)

num_sols = length(w)/length(PUApprox);

g_f  = JacobianFowardLU(PUApprox,g_mat.l,g_mat.u,g_mat.p,w);
 
g_f_hat = [];
for i=1:num_sols
    g_f_hat = [g_f_hat;PUApprox.Fine2Coarse(g_f((1:length(PUApprox)) + (i-1)*length(PUApprox)),j)];
end

w_hat = [];
for i=1:num_sols
    w_hat = [w_hat;PUApprox.Fine2Coarse(w((1:length(PUApprox)) + (i-1)*length(PUApprox)),j)];
end

PUApprox.Coarsen();

g_hat = JacobianFowardLU(PUApprox,g_p_mat.l,g_p_mat.u,g_p_mat.p,w_hat); 

R = g_hat - g_f_hat;

b = -(J_er+T_hat)*R;

er = J_er\b;

er = er - w_hat;

r_er = [];

for i=1:num_sols
    r_er = [r_er;PUApprox.Coarse2Fine(er((1:length(PUApprox)) + (i-1)*length(PUApprox)))];
end

PUApprox.Refine();

end


