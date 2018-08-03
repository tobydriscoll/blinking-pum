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
function [z,l,u,p,J_v_pls_er] = ParPreconditionedTwoLevel(sol,PUApprox,evalF,Jac)

[c_sol,J_v_pls_er ] = CoarseCorrect(PUApprox,sol,evalF,Jac);

[z,l,u,p] = ParPreconditionedNewtonForward(sol+c_sol,PUApprox,evalF,Jac);

z  = c_sol + z;

end





