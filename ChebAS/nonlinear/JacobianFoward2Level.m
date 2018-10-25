% JacobianFoward2Level
% Matrix free evaluation for the Jacobian of the SNK method.
%
%
% INPUT:
%      PUApprox: PUApprox approximation        
%        L,U,p: cellarray of LU matrices and permutation vector for
%               patches of f(v+e)
%      Jac_hat: Jacobian used for the nonlinear coarse correction
%          FJv: cell array of local jacobians on fine grid
%       FJ_hat: cell array of local jacobians on the coarse grid
%            w: solution
%            j: difference by power of two between coarse and fine grids.
%               Grids are chosen from 2,5,9,33... . For example if we
%               have coarse and fine grids of 9,33 (in each dimension)
%               then j=2. (This should just be apart of PUApprox).
%
% NOTE x is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then x = [u1;u2;v1;v2].
function [ y ] = JacobianFoward2Level(PUApprox,L,U,p,Jac_hat,FJv,FJv_hat,w,j)

c_w = LinearCoarseCorrect( PUApprox, w,Jac_hat,FJv,FJv_hat,j);

y = c_w + JacobianFowardLU(PUApprox,L,U,p,w+c_w);

end

