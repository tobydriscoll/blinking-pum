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
function [z,l,u,p,J_v_pls_er] = ParPreconditionedTwoLevel3(sol,PUApprox,Leaf,G,evalF,Jac,tol,tol_c)

[c_sol,J_v_pls_er ] = CoarseCorrect3(PUApprox,Leaf,G,sol,evalF,Jac,tol_c);

[z,l,u,p] = ParPreconditionedNewtonForward(sol+c_sol,PUApprox,evalF,Jac,tol);

z  = c_sol + z;

end





