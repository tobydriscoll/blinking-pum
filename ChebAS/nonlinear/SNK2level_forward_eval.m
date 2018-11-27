% SNK_forward_eval
%
% The residual minimized for The Schwarz Newton Krylov (SNK) method.
%
% INPUT:
%
%            sol: solution used for computing the residual.
%
%            PUApprox: PUApprox approximation.
%
%            evalF: nonlinear residual function f(x,p) for solution x and local
%             approximation p. f(x,p) evaluates the residual on the domain
%             of p.
%
%            Jac: Jacobian function f(x,p) for solution x and local
%             approximation p. f(x,p) evaluates the residual on the domain
%             of p.
%
%              j: difference by power of two between coarse and fine grids.
%                 Grids are chosen from 2,5,9,33... . For example if we
%                 have coarse and fine grids of 9,33 (in each dimension)
%                 then j=2. (This should just be apart of PUApprox).
%
%          tol_n: [rel_tol, abs_tol] relative and absolute tolerance used
%                 for local Newtons method.
%
%          tol_c: tolerance used for nonlinear coarse correction.
%
% OUTPUT:
%          z: correction of solution
%
%      l,u,p: cell array of LU decomposition and ordering for local
%             Jacobians.
%
% J_v_pls_er: Local sparse Jacobian used in solving the nonlinear coarse
%             correction
%
% NOTE sol is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then sol = [u1;u2;v1;v2].
function [z,l,u,p,J_v_pls_er,c_sol] = SNK2level_forward_eval(sol,PUApprox,evalF,Jac,j,tol,tol_c)

[c_sol,J_v_pls_er ] = CoarseCorrect(PUApprox,sol,evalF,Jac,j,tol_c);

[z,l,u,p] = SNK_forward_eval(sol+c_sol,PUApprox,evalF,Jac,tol);

z  = c_sol + z;
end





