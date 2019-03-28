% JacobianFowardLU
% Matrix free evaluation for the Jacobian of the SNK method.
%
% INPUT:
%    PUApproxArray: cell array of PUApprox approximations
%            L,U,p: cell array of the LU decomposition of the local Jacobians.
%                x: solution
%
% NOTE x is presumed to be ordered by solution first, then patch.
%      For example, suppose there are two patches p1, p2 each with
%      two solutions u1 v1, u2 v2. Then x = [u1;u2;v1;v2].
%
% ALSO NOTE This method takes into account of boundary info is accounted
% for by the effect if the individual trees are packed or not.
function [ output ] = BlockLinearResidual(PUApproxArray,J,x)

if ~iscell(PUApproxArray)
    PUApproxArray = {PUApproxArray};
end

num_leaves = length(PUApproxArray{1}.leafArray);

x = unpackPUvecs(x,PUApproxArray);

output = x;

%Parallel part
for k=1:num_leaves
    
    output{k} =  J{k}*x{k};
    
end

%Take {[u1;u2],[v1;v2]} to [u1;u2;v1;v2]
output = packPUvecs(output,PUApproxArray);
