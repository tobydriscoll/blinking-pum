% INPUT:
%      PUApprox: PUApprox approximation
%           sol: given solution
%         evalF: residual function for leaf
%           Jac: jacobian function for leaf
%
% OUTPUT:
%          z: correction of solution
%          J: cell array of local Jacobians
%
% NOTE sol is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then sol = [u1;u2;v1;v2].
function [z,l,u,p,g_mat,g_p_mat,J_er] = ParPreconditionedTwoLevelG(sol,PUApprox,evalF,Jac,tol,tol_c,j)

[ er, g_mat,g_p_mat,J_er] = CoarseCorrect2(PUApprox,sol,evalF,Jac,tol,tol_c,j);

[z,l,u,p] = ParPreconditionedNewtonForward(sol+er,PUApprox,evalF,Jac,tol);

z  = er + z;

end
