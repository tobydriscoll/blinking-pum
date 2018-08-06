% INPUT:
%      PUApprox: PUApprox approximation        
%        L,U,p: cellarray of LU matrices and permutation vector for
%               patches of f(v+e)
%      Jac_hat: f_hat'(v+e)-T_hat
%        T_hat: T_hat
%            v: solution v
%          rhs: rhs to be used with gmres
%
% NOTE x is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then x = [u1;u2;v1;v2].
function [ y ] = GMRES_2_level(PUApprox,jac_f,L,U,p,Jac_hat,v,rhs)

num_sols = length(v)/length(PUApprox);

T_hat  = CoarseInterfaceInterp(PUApprox,num_sols);

[FJv,FJv_hat] = ComputeJac(PUApprox,jac_f,v);

FJv_hat = blkdiag(FJv_hat{:});

[y,~,~,~,~] = gmres(@(w)JacobianFoward2Level(PUApprox,L,U,p,Jac_hat,T_hat,FJv,FJv_hat,w),rhs,[],1e-14);
end