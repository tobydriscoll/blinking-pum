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
function [ output ] = ASPreconditionerTime(PUApproxArray,L,U,p,x,sub_ind,loc_sub_ind)

if ~iscell(PUApproxArray)
    PUApproxArray = {PUApproxArray};
end

take_sub_ind = nargin>6;

if take_sub_ind
    x_temp = zeros(length(sub_ind),1);
    x_temp(sub_ind) = x;
    x = x_temp;
end

num_leaves = length(PUApproxArray{1}.leafArray);
x = unpackPUvecs(x,PUApproxArray);

%Parallel part
for k=1:num_leaves
    if take_sub_ind
        x_t = x{k}(loc_sub_ind{k});
        output{k} = zeros(size(x{k}));
        output{k}(loc_sub_ind{k}) = U{k}\(L{k}\x_t(p{k}));
    else
        output{k} = U{k}\(L{k}\x{k}(p{k}));
    end
end

%Take {[u1;u2],[v1;v2]} to [u1;u2;v1;v2]
output = packPUvecs(output,PUApproxArray);

if take_sub_ind
    output = output(sub_ind);
end