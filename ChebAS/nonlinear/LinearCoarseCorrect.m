% LinearCoarseCorrect
% This methods computes the linear coarse correction for the two level SNK
% method.
%
% INPUT:
%      PUApprox: PUApprox approximation        
%            w: the solution.
%
% Lc, Uc Pc Qc: LU decomposition with row, column permutation matrices of
%                jacoabian used for coarse solve
%
%          FJv: cell array of local jacobians on fine grid.
%
%       FJ_hat: cell array of local jacobians on the coarse grid.
%
%            j: difference by power of two between coarse and fine grids.
%               Grids are chosen from 2,5,9,33... . For example if we
%               have coarse and fine grids of 9,33 (in each dimension)
%               then j=2. (This should just be apart of PUApprox).
%
% OUTPUT:
%            y: linear coarse correction.
%
% NOTE x is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then x = [u1;u2;v1;v2].
function [ y ] = LinearCoarseCorrect( PUApprox, w,Lc,Uc,Pc,Qc,FJv,FJv_hat,j)


num_sols = length(w)/length(PUApprox);

w_hat = [];
for i=1:num_sols
    w_hat = [w_hat;PUApprox.Fine2Coarse(w((1:length(PUApprox)) + (i-1)*length(PUApprox)),j)];
end

%need to linearize this
r = ParLinearResidual(w,PUApprox,FJv);

r_hat = [];

for i=1:num_sols
    r_hat = [r_hat;PUApprox.Fine2Coarse(r((1:length(PUApprox)) + (i-1)*length(PUApprox)),j)];
end

PUApprox.Coarsen();

b_hat = ParLinearResidual(w_hat,PUApprox,FJv_hat)-r_hat;

y_hat = Qc*(Uc\(Lc\(Pc*b_hat)))-w_hat;

y = [];

for i=1:num_sols
    y = [y;PUApprox.Coarse2Fine(y_hat((1:length(PUApprox)) + (i-1)*length(PUApprox)))];
end

PUApprox.Refine();

end
