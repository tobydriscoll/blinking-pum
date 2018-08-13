%      PUApprox: PUApprox approximation        
%        L,U,p: cellarray of LU matrices and permutation vector for
%               patches of f(v+e)
%      Jac_hat: f_hat'(v+e)-T_hat
%        T_hat: T_hat
%          FJv: cell array of f'(v_j)
%       FJ_hat: sparse matrix of f_hat'(v_hat)
%            w: solution
%
% NOTE x is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then x = [u1;u2;v1;v2].
function [ y ] = JacobianFoward2LevelG(PUApprox,L,U,p,g_mat,g_p_mat,J_er,T_hat,w,j)

c_w = LinearCoarseCorrect2(PUApprox,w,g_mat,g_p_mat,J_er,T_hat,j);

y = c_w + JacobianFowardLU(PUApprox,L,U,p,w+c_w);

end