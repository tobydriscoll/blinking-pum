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
function [ output ] = JacobianForwardLUTime(PUApproxArray,L,U,p,x,alpha)

if ~iscell(PUApproxArray)
    PUApproxArray = {PUApproxArray};
end

num_sols = length(PUApproxArray);

num_leaves = length(PUApproxArray{1}.leafArray);

sol_lengths = zeros(num_sols,1);

for i=1:num_sols
    sol_lengths(i) = length(PUApproxArray{i});
end

x = mat2cell(x,sol_lengths);

x_unpacked = x;

for i=1:num_sols
    
    if PUApproxArray{i}.is_packed
        
        % pull the boundry info for the packed functions. The boundary
        % info is stored within the tree itself. 
        % Set boundary to zero (This is important!)
        x_unpacked{i} = PUApproxArray{i}.Getunpackedvalues(x{i},true);
        
    else
        
        x_unpacked{i} = x{i};
        
    end
    
end

%Take [u1;u2;v1;v2] to {[u1;u2],[v1;v2]}
x = cell2mat(x);
[x,lens] = unpackPUvecs(x,PUApproxArray);

z = x;

%non parallel part
for k=1:num_leaves
    
    z_loc = mat2cell(zeros(size(x{k})),lens{k});
    
    for i=1:num_sols
        in_border = PUApproxArray{i}.leafArray{k}.inner_boundary;
        z_loc{i}(in_border) = alpha*PUApproxArray{i}.leafArray{k}.Binterp*x_unpacked{i};
    end
    
    z{k} = cell2mat(z_loc);
    
end

output = x;

%Parallel part
for k=1:num_leaves
    
    output{k} = U{k}\(L{k}\z{k}(p{k})) - x{k};
    
   % output{k} =  L{k}*(U{k}*z{k}(p{k}))*x{k}-z{k};
    
end

%Take {[u1;u2],[v1;v2]} to [u1;u2;v1;v2]
output = packPUvecs(output,PUApproxArray);